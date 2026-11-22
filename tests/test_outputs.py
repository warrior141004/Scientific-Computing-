import hashlib
import os

import grader


PUBLIC_HASHES = {
    "calibration.tsv": "2350a5210c6b7cf6db3ef3e182b3bfa8a88f7b4c969c79821b39f587b93b2278",
    "clock_beacons.tsv": "0da5b738515dac29c40b719d5a73a1ad44dcfb6c20a7422bcdf6dbf4cc382f49",
    "constants.cfg": "d82876c17458928719a48d31eba755b1e5a5ed46706afc09371c67467941827a",
    "manifest.txt": "be9109299fc42f1b2562e477b18634e209665077de2fa32fe143203640573840",
    "model.cfg": "13cb40645ce1e8f90ebab60294d41de0b6506b84f47201416a81407ea45ee5b3",
    "observations.tsv": "e2c94100f1456e56edc74527f9f2c413db727b1276cef6209f3209fd4df12fa7",
    "queries.tsv": "d3e525ba3cc774629a36299a58b2114e66b5be6dbc1d81bf5a6231a18bd4713f",
    "stations.tsv": "57d36197db6fec42998ff7090aac6cbd7e1deb31cb375a83da566906a90bca33",
}


def test_public_input_integrity():
    """The public fixture in /app/input_data/public should be unchanged."""
    public_dir = os.path.join(grader.APP_DIR, "input_data", "public")
    for filename, expected in PUBLIC_HASHES.items():
        with open(os.path.join(public_dir, filename), "rb") as f:
            actual = hashlib.sha256(f.read()).hexdigest()
        assert actual == expected


def test_all_bundles_pass():
    """The solver should pass public and compatible hidden-style bundles."""
    reports = {}
    for name, bundle, truth in grader.bundles():
        reports[name] = grader.evaluate(bundle, truth)
        assert reports[name]["ok"], "%s: %s" % (name, reports[name])


def test_generated_matrix_is_present():
    """The verifier should include deterministic generated interaction cases."""
    generated = [(name, bundle, truth) for name, bundle, truth in grader.bundles() if name.startswith("generated_")]
    assert len(generated) >= 10
    cubic = 0
    zero_response = 0
    for _name, bundle, truth in generated:
        inputs = grader.manifest_paths(bundle)
        constants = {row[0]: row[1] for row in grader.rows(inputs["constants"]) if len(row) >= 2}
        if constants.get("clock_order") == "3":
            cubic += 1
        for row in grader.rows(os.path.join(truth, "truth_gains.tsv")):
            if len(row) >= 5 and abs(float(row[4])) > 0.0:
                zero_response += 1
    assert cubic >= 6
    assert zero_response >= len(generated)


def test_hidden_bundles_are_distinct():
    """The included holdout bundles should not be duplicate copies."""
    a = open(grader.bundles()[1][1] + "/stations.tsv", encoding="utf-8").read()
    b = open(grader.bundles()[2][1] + "/stations.tsv", encoding="utf-8").read()
    assert a != b
