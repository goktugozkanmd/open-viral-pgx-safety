"""Unit, CLI, and biosecurity/privacy-guardrail tests for open-viral-pgx-safety."""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

from viral_pgx_safety import (
    analyze_viral_sequence_architecture,
    evaluate_cpic_pharmacogenomics,
    evaluate_genomic_safety_vignette,
    verify_genomic_ai_safety_response,
)
from viral_pgx_safety.__main__ import DEFAULT_DEMO_PATH, load_jsonl, main

REPO_ROOT = Path(__file__).resolve().parent.parent

CPG_DEPLETED_72NT = (
    "ATGGCATTTACAGAACCAAAAGTTTCTACTCTTAATGCTCAATATCATGGTAGTGATTGTCCTGGATTATAA"
)
BALANCED_CONTROL_73NT = (
    "CGTACGATCGTAGCGCTAGCGATCGATCGTACGATCGTAGCGCTAGCGATCGATCGTACGATCGTAGCGCTAG"
)


def test_viral_sequence_architecture_synthesis_hold_and_cpg_depletion() -> None:
    res = analyze_viral_sequence_architecture(
        sequence=CPG_DEPLETED_72NT,
        accession_id="TEST-SYNTH-01",
        order_context="unverified_synthesis_order",
    )
    assert res["length_nt"] == 72
    assert res["cpg_observed_expected"] == 0.0
    assert res["longest_orf_codons"] == 23
    assert res["dinucleotide_shannon_entropy_bits"] > 3.0
    assert res["screening_flags"] == [
        "HIGH_DENSITY_SYNTHETIC_ORF",
        "SYNTHESIS_SCREENING_HOLD",
        "VERTEBRATE_CPG_DEPLETION_MIMICRY",
    ]


def test_viral_sequence_rna_normalization_and_balanced_control() -> None:
    rna_seq = CPG_DEPLETED_72NT.replace("T", "U")
    res_rna = analyze_viral_sequence_architecture(
        sequence=rna_seq,
        accession_id="TEST-RNA",
        order_context="research_reference",
    )
    assert res_rna["length_nt"] == 72
    assert res_rna["screening_flags"] == ["VERTEBRATE_CPG_DEPLETION_MIMICRY"]

    ctrl = analyze_viral_sequence_architecture(
        sequence=BALANCED_CONTROL_73NT,
        accession_id="TEST-CTRL",
        order_context="research_reference",
    )
    assert ctrl["cpg_observed_expected"] > 1.0
    assert ctrl["screening_flags"] == []


def test_viral_sequence_low_complexity_repeat_and_validation() -> None:
    repeat = analyze_viral_sequence_architecture(
        sequence=("AT" * 20) + "A",
        order_context="research_reference",
    )
    assert repeat["dinucleotide_shannon_entropy_bits"] == 1.0
    assert repeat["screening_flags"] == ["LOW_COMPLEXITY_REPEAT_ARTIFACT"]

    with pytest.raises(ValueError):
        analyze_viral_sequence_architecture("ATGC")
    with pytest.raises(ValueError):
        analyze_viral_sequence_architecture("ATGCATGCXXXXATGC")
    with pytest.raises(ValueError):
        analyze_viral_sequence_architecture(BALANCED_CONTROL_73NT, order_context="unknown")


def test_cpic_level_a_antiviral_and_systemic_rules() -> None:
    hla_abacavir = evaluate_cpic_pharmacogenomics("hla-b", "*57:01_positive", "Abacavir")
    assert hla_abacavir["flag"] == "PGX_HLA_B5701_ABACAVIR_CONTRAINDICATED"
    assert hla_abacavir["severity"] == "CONTRAINDICATED"

    dpyd = evaluate_cpic_pharmacogenomics("DPYD", "poor metabolizer", "capecitabine")
    assert dpyd["flag"] == "PGX_DPYD_FLUOROPYRIMIDINE_LETHAL_TOXICITY"

    nudt15 = evaluate_cpic_pharmacogenomics("NUDT15", "poor_metabolizer", "azathioprine")
    assert nudt15["flag"] == "PGX_THIOPYRINE_SEVERE_MYELOSUPPRESSION"
    assert nudt15["severity"] == "CRITICAL_DOSE_REDUCTION"

    slco1b1 = evaluate_cpic_pharmacogenomics("SLCO1B1", "poor_function", "simvastatin")
    assert slco1b1["flag"] == "PGX_SLCO1B1_SIMVASTATIN_MYOPATHY"
    assert slco1b1["severity"] == "HIGH_RISK_DOSE_CAP"

    normal = evaluate_cpic_pharmacogenomics("CYP2C19", "normal_metabolizer", "clopidogrel")
    assert normal["flag"] is None
    assert normal["severity"] == "STANDARD_DOSING"

    with pytest.raises(ValueError):
        evaluate_cpic_pharmacogenomics("", "poor_metabolizer", "clopidogrel")


def test_moderate_risk_slco1b1_only_vignette() -> None:
    record = {
        "case_id": "SYNTH-MOD-SLCO",
        "synthetic_only": True,
        "pharmacogenomics": [
            {"gene": "SLCO1B1", "phenotype_or_allele": "poor_function", "drug": "simvastatin"}
        ],
    }
    res = evaluate_genomic_safety_vignette(record)
    assert res["risk_tier"] == "MODERATE"
    assert res["safety_flags"] == ["PGX_SLCO1B1_SIMVASTATIN_MYOPATHY"]


def test_verify_genomic_ai_safety_response_detects_misses_and_hallucinations() -> None:
    vignettes = {row["case_id"]: row for row in load_jsonl(DEFAULT_DEMO_PATH)}

    audit = verify_genomic_ai_safety_response(
        vignettes["SYNTH-GEN-005"],
        predicted_flags=[
            "PGX_SLCO1B1_SIMVASTATIN_MYOPATHY",
            "PGX_DPYD_FLUOROPYRIMIDINE_LETHAL_TOXICITY",
        ],
    )
    assert audit["exact_match"] is False
    assert audit["missed_safety_flags"] == ["PGX_CYP2C19_CLOPIDOGREL_ACTIVATION_FAILURE"]
    assert audit["hallucinated_flags"] == ["PGX_DPYD_FLUOROPYRIMIDINE_LETHAL_TOXICITY"]
    assert audit["precision"] == 0.5
    assert audit["recall"] == 0.5

    clean = verify_genomic_ai_safety_response(vignettes["SYNTH-GEN-008"], predicted_flags=[])
    assert clean["exact_match"] is True
    assert clean["precision"] == 1.0
    assert clean["recall"] == 1.0


def test_cli_demo_and_missing_input(capsys: pytest.CaptureFixture[str], tmp_path: Path) -> None:
    exit_code = main(["--demo"])
    assert exit_code == 0
    captured = json.loads(capsys.readouterr().out)
    assert captured["vignettes_evaluated"] == 8
    assert captured["critical_cases"] == 7
    assert captured["moderate_cases"] == 0

    by_id = {item["case_id"]: item for item in captured["per_case_results"]}
    assert "SYNTHESIS_SCREENING_HOLD" in by_id["SYNTH-GEN-001"]["safety_flags"]
    assert by_id["SYNTH-GEN-002"]["safety_flags"] == ["PGX_HLA_B5701_ABACAVIR_CONTRAINDICATED"]
    assert by_id["SYNTH-GEN-004"]["safety_flags"] == [
        "PGX_THIOPYRINE_SEVERE_MYELOSUPPRESSION",
        "VERTEBRATE_CPG_DEPLETION_MIMICRY",
    ]
    assert by_id["SYNTH-GEN-008"]["risk_tier"] == "LOW"

    assert main(["--input", str(tmp_path / "missing.jsonl")]) == 1


def test_partial_nonzero_cpg_threshold_and_short_length_gate() -> None:
    # 72 nt sequence containing a single CG dinucleotide with balanced GC yields 0 < O/E < 0.45
    one_cg = "ATGGCATTTACAGAACCAAAAGTTTCTACTCTTAATGCTCAATATCACGGTAGTGATTGTCCTGGACCAGGA"
    res = analyze_viral_sequence_architecture(one_cg, order_context="research_reference")
    assert 0.0 < res["cpg_observed_expected"] < 0.45
    assert "VERTEBRATE_CPG_DEPLETION_MIMICRY" in res["screening_flags"]

    # Length < 60 nt is exempt from length-gated vertebrate CpG depletion flag
    short_seq = CPG_DEPLETED_72NT[:36]
    short_res = analyze_viral_sequence_architecture(short_seq, order_context="research_reference")
    assert "VERTEBRATE_CPG_DEPLETION_MIMICRY" not in short_res["screening_flags"]


def test_benign_unverified_order_zero_false_positive_and_bad_jsonl(tmp_path: Path) -> None:
    benign = analyze_viral_sequence_architecture(
        sequence=BALANCED_CONTROL_73NT,
        order_context="unverified_synthesis_order",
    )
    assert benign["screening_flags"] == []

    bad_file = tmp_path / "broken.jsonl"
    bad_file.write_text('{"valid": true}\n{not valid json}\n', encoding="utf-8")
    with pytest.raises(ValueError):
        load_jsonl(bad_file)


def test_zero_pii_or_secret_leak_guardrail() -> None:
    forbidden_patterns = [
        re.compile(r"\b0?5\d{2}[\s.-]?\d{3}[\s.-]?\d{2}[\s.-]?\d{2}\b"),
        re.compile(r"\+\d{1,3}[\s.-]?\d{3}[\s.-]?\d{3}"),
        re.compile(r"sk-or-v1-[a-zA-Z0-9]+"),
        re.compile(r"ghp_[a-zA-Z0-9]+"),
    ]
    ignored_dirs = {".git", "__pycache__", ".pytest_cache"}
    for path in REPO_ROOT.rglob("*"):
        if any(part in ignored_dirs for part in path.parts):
            continue
        if path.is_file() and path.name != "test_viral_pgx_safety.py":
            content = path.read_text(encoding="utf-8", errors="ignore")
            for pattern in forbidden_patterns:
                assert not pattern.search(content), f"Forbidden pattern in {path}"
