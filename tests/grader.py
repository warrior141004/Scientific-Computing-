import atexit
import cmath
import math
import os
import shutil
import subprocess
import tempfile

APP_DIR = os.environ.get("APP_DIR", "/app")
HERE = os.path.dirname(os.path.abspath(__file__))
MU_TOL = 2.0e-6
PHASE_TOL = 0.30
LOG_TOL = 0.0030
GAIN_TOL = 2.0e-7
RUN_TIMEOUT = 120
COMPILE_TIMEOUT = 120
_COMPILED_DIR = None
_COMPILED_EXE = None
_GENERATED_ROOT = None


def cdiff(a, b):
    return (a - b + 180.0) % 360.0 - 180.0


def frequency_from_tick(tick, clock, beacon_bounds=None):
    def value(x):
        total = 0.0
        power = 1.0
        for coef in clock:
            total += coef * power
            power *= x
        return total

    def slope(x):
        total = 0.0
        power = 1.0
        for i, coef in enumerate(clock[1:], start=1):
            total += i * coef * power
            power *= x
        return total

    roots = []
    lo, hi = -20.0, 20.0
    scans = 4000
    px = lo
    py = value(px) - tick
    for i in range(1, scans + 1):
        x = lo + (hi - lo) * i / scans
        y = value(x) - tick
        if py == 0.0 or y == 0.0 or (py < 0.0) != (y < 0.0):
            a, b, fa = px, x, py
            for _ in range(90):
                m = 0.5 * (a + b)
                fm = value(m) - tick
                if (fa < 0.0) == (fm < 0.0):
                    a, fa = m, fm
                else:
                    b = m
            roots.append(0.5 * (a + b))
        px, py = x, y
    valid = [x for x in roots if slope(x) > 0.0 and -20.0 < x < 20.0]
    if beacon_bounds is None:
        chosen = min(valid, key=abs)
    else:
        blo, bhi = beacon_bounds
        center = 0.5 * (blo + bhi)
        chosen = min(valid, key=lambda x: ((10.0 if x < blo - 1e-9 or x > bhi + 1e-9 else 0.0)
                                           + abs(x - center) * 1e-3
                                           + abs(sum(coef * (x ** i) for i, coef in enumerate(clock)) - tick)))
    return math.exp(chosen)


def rows(path):
    out = []
    with open(path, encoding="utf-8") as f:
        for raw in f:
            s = raw.strip()
            if not s or s.startswith("#"):
                continue
            out.append(s.split())
    return out


def manifest_paths(bundle):
    out = {}
    for r in rows(os.path.join(bundle, "manifest.txt")):
        if len(r) >= 2:
            out[r[0]] = os.path.join(bundle, r[1])
    return out


def copy_bundle(bundle, sandbox):
    target = os.path.join(sandbox, "bundle")
    shutil.copytree(bundle, target)
    return target


def compiled_solver():
    global _COMPILED_DIR, _COMPILED_EXE
    if _COMPILED_EXE is not None:
        return _COMPILED_EXE
    src = os.path.join(APP_DIR, "solve.cpp")
    if not os.path.exists(src):
        raise RuntimeError("expected solver source at /app/solve.cpp")
    _COMPILED_DIR = tempfile.mkdtemp(prefix="colecole_build_")
    _COMPILED_EXE = os.path.join(_COMPILED_DIR, "solve")
    atexit.register(lambda: shutil.rmtree(_COMPILED_DIR, ignore_errors=True))
    p = subprocess.run(["g++", "-std=c++17", "-O2", "-pipe", src, "-o", _COMPILED_EXE],
                       capture_output=True, text=True, timeout=COMPILE_TIMEOUT)
    if p.returncode:
        raise RuntimeError("compile failed: " + p.stderr[-2000:])
    return _COMPILED_EXE


def run_program(bundle):
    sandbox = tempfile.mkdtemp(prefix="colecole_run_")
    try:
        copied = copy_bundle(bundle, sandbox)
        p = subprocess.run([compiled_solver(), os.path.join(copied, "manifest.txt")], cwd=sandbox,
                           capture_output=True, text=True, timeout=RUN_TIMEOUT)
        if p.returncode:
            raise RuntimeError("solver exited %d: %s" % (p.returncode, p.stderr[-2000:]))
        return sandbox
    except Exception:
        shutil.rmtree(sandbox, ignore_errors=True)
        raise


def parse_truth(path):
    model = {}
    for r in rows(os.path.join(path, "truth_model.tsv")):
        model[r[0]] = float(r[1])
    gains = {
        r[0]: (float(r[1]), float(r[2]), float(r[3]), float(r[4]) if len(r) > 4 else 0.0)
        for r in rows(os.path.join(path, "truth_gains.tsv"))
    }
    clock = tuple(float(x) for x in rows(os.path.join(path, "truth_clock.tsv"))[0])
    pred = {int(r[0]): (float(r[1]), float(r[2])) for r in rows(os.path.join(path, "truth_predictions.tsv"))}
    query = {int(r[0]): (float(r[1]), float(r[2])) for r in rows(os.path.join(path, "truth_queries.tsv"))}
    return model, gains, clock, pred, query


def wrap_phase(deg):
    return (deg + 180.0) % 360.0 - 180.0


def forward_model(freq, model, constants):
    pi = constants["pi"]
    omega_i = 1j * 2.0 * pi * freq
    sigma = [1.0 / model["rho0"], 1.0 / model["rho1"], 1.0 / model["rho2"]]
    d = int(model["disp"])
    cole = (omega_i * model["tau"]) ** model["c"]
    sigma[d] *= 1.0 + model["m"] / (1.0 + cole)
    z = cmath.sqrt(omega_i * constants["mu0"] / sigma[2])
    for j, hkey in ((1, "h1"), (0, "h0")):
        gamma = cmath.sqrt(omega_i * constants["mu0"] * sigma[j])
        eta = cmath.sqrt(omega_i * constants["mu0"] / sigma[j])
        t = cmath.tanh(gamma * model[hkey])
        z = eta * (z + eta * t) / (eta + z * t)
    return z


def log_polar(z, pi):
    return math.log10(abs(z)), math.atan2(z.imag, z.real) * 180.0 / pi


def write_lines(path, lines):
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
        f.write("\n")


def tick_from_freq(freq, clock):
    x = math.log(freq)
    total = 0.0
    power = 1.0
    for coef in clock:
        total += coef * power
        power *= x
    return total


def generated_case_specs():
    return [
        (0, (70, 380, 1600), (38, 125), 0.29, 0.14, 0.68, (8.7, 2.20, 0.055), False),
        (1, (95, 520, 2100), (52, 160), 0.36, 0.22, 0.74, (9.4, 2.05, -0.038, 0.0060), False),
        (2, (55, 280, 1450), (31, 95), 0.24, 0.09, 0.61, (7.8, 2.35, 0.072, -0.0045), True),
        (0, (130, 760, 2500), (70, 190), 0.41, 0.31, 0.82, (10.2, 1.92, -0.026, 0.0038), False),
        (1, (45, 610, 1700), (28, 210), 0.33, 0.055, 0.58, (8.1, 2.48, 0.044), True),
        (2, (160, 430, 3200), (88, 145), 0.27, 0.47, 0.77, (11.0, 1.85, -0.041, 0.0055), False),
        (0, (82, 260, 1200), (46, 80), 0.48, 0.18, 0.70, (6.9, 2.65, 0.030, -0.0032), False),
        (1, (210, 920, 3600), (115, 260), 0.22, 0.72, 0.86, (12.5, 1.70, -0.018), True),
        (2, (62, 700, 1850), (24, 175), 0.44, 0.12, 0.64, (9.0, 2.16, 0.063, 0.0042), False),
        (0, (115, 470, 2800), (64, 118), 0.31, 0.39, 0.79, (10.8, 1.98, -0.032, -0.0028), True),
    ]


def generate_bundle(root, idx, spec):
    disp, rho, thick, m, tau, cole_c, clock, nested = spec
    label = "generated_%02d" % idx
    bundle = os.path.join(root, label, "bundle")
    truth = os.path.join(root, label, "truth")
    data_dir = os.path.join(bundle, "data") if nested else bundle
    os.makedirs(data_dir, exist_ok=True)
    os.makedirs(truth, exist_ok=True)

    constants = {"mu0": 4.0 * math.pi * 1e-7, "pi": math.pi, "huber_delta": 0.018}
    model = {
        "rho0": float(rho[0]), "rho1": float(rho[1]), "rho2": float(rho[2]),
        "h0": float(thick[0]), "h1": float(thick[1]),
        "m": m, "tau": tau, "c": cole_c, "disp": disp,
    }
    stations = [
        ("S%02dA" % idx, -172.0 + 9.0 * idx, complex(1.04 + 0.012 * idx, 0.10 - 0.015 * idx), 0.0020 + 0.00018 * (idx % 4), 0.00035 + 0.00004 * (idx % 5)),
        ("S%02dB" % idx, 37.0 + 11.0 * idx, complex(0.88 + 0.020 * idx, -0.08 + 0.010 * idx), 0.0042 + 0.00015 * (idx % 5), 0.00085 + 0.00005 * (idx % 4)),
        ("S%02dC" % idx, 196.0 - 7.0 * idx, complex(1.18 - 0.010 * idx, 0.04 + 0.012 * idx), 0.0016 + 0.00022 * (idx % 3), 0.00015 + 0.000035 * (idx % 6)),
        ("S%02dD" % idx, -41.0 + 5.0 * idx, complex(0.97 + 0.006 * idx, 0.16 - 0.009 * idx), 0.0031 + 0.00012 * (idx % 6), 0.00065 + 0.00005 * (idx % 3)),
    ]

    prefix = "data/" if nested else ""
    write_lines(os.path.join(bundle, "manifest.txt"), [
        "constants %sconstants.cfg" % prefix,
        "model %smodel.cfg" % prefix,
        "clock %sclock_beacons.tsv" % prefix,
        "stations %sstations.tsv" % prefix,
        "calibration %scalibration.tsv" % prefix,
        "observations %sobservations.tsv" % prefix,
        "queries %squeries.tsv" % prefix,
    ])
    write_lines(os.path.join(data_dir, "constants.cfg"), [
        "mu0 %.17g" % constants["mu0"],
        "pi %.17g" % constants["pi"],
        "huber_delta %.17g" % constants["huber_delta"],
        "clock_order %d" % (len(clock) - 1),
    ])
    write_lines(os.path.join(data_dir, "model.cfg"), [
        "layer_count 3",
        "dispersive_layer %d" % disp,
        "initial_rho %.17g %.17g %.17g" % (rho[0] * 0.83, rho[1] * 1.19, rho[2] * 0.91),
        "initial_thickness %.17g %.17g" % (thick[0] * 1.22, thick[1] * 0.84),
        "initial_m %.17g" % min(0.78, max(0.02, m * 0.82 + 0.04)),
        "initial_tau %.17g" % max(0.003, tau * 1.35),
        "initial_c %.17g" % min(0.98, max(0.22, cole_c * 0.90)),
        "lower_rho 10 10 50",
        "upper_rho 5000 5000 8000",
        "lower_thickness 5 5",
        "upper_thickness 500 700",
        "lower_m 0.01",
        "upper_m 0.80",
        "lower_tau 0.002",
        "upper_tau 5.0",
        "lower_c 0.20",
        "upper_c 1.00",
    ])
    write_lines(os.path.join(data_dir, "stations.tsv"),
                ["# station_id orientation_degrees"] + ["%s\t%.17g" % (s[0], s[1]) for s in stations])
    beacon_freqs = [0.045, 0.09, 0.18, 0.42, 1.1, 2.9, 7.6, 19.0, 48.0]
    write_lines(os.path.join(data_dir, "clock_beacons.tsv"),
                ["# beacon_id instrument_tick known_frequency_hz"] +
                ["B%d\t%.17g\t%.17g" % (i, tick_from_freq(f, clock), f) for i, f in enumerate(beacon_freqs)])

    cal_freqs = [0.055, 0.14, 0.37, 0.95, 2.6, 6.8, 18.0, 46.0, 92.0]
    cal_lines = ["# station_id frequency reference_re reference_im measured_re measured_im"]
    for si, (sid, _angle, gain, pole, zero) in enumerate(stations):
        for j, freq in enumerate(cal_freqs):
            ref = complex(0.94 + 0.035 * j + 0.012 * si, -0.24 + 0.045 * j - 0.006 * si)
            measured = gain * complex(1.0, 2.0 * math.pi * freq * zero) / complex(1.0, 2.0 * math.pi * freq * pole) * ref
            cal_lines.append("%s\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g" %
                             (sid, freq, ref.real, ref.imag, measured.real, measured.imag))
    write_lines(os.path.join(data_dir, "calibration.tsv"), cal_lines)

    obs_freqs = [math.exp(-2.45 + 0.285 * j) for j in range(23)]
    obs_lines = ["# observation_id station_id instrument_tick raw_re raw_im"]
    truth_pred_lines = []
    oid_base = 1000 + idx * 200
    outlier_index = (idx * 7 + 5) % (len(stations) * len(obs_freqs))
    obs_index = 0
    for si, (sid, angle, gain, pole, zero) in enumerate(stations):
        orient = cmath.exp(1j * angle * math.pi / 180.0)
        for j, freq in enumerate(obs_freqs):
            z = forward_model(freq, model, constants)
            la, ph = log_polar(z, constants["pi"])
            wobble = 0.00055 * complex(math.sin((idx + 1) * (j + 2) * 0.73),
                                       math.cos((si + 2) * (j + 1) * 0.61))
            z_obs = z * (1.0 + wobble)
            if obs_index == outlier_index:
                z_obs *= complex(1.16 - 0.015 * idx, -0.13 + 0.02 * (idx % 3))
            raw = z_obs * orient * gain * complex(1.0, 2.0 * math.pi * freq * zero) / complex(1.0, 2.0 * math.pi * freq * pole)
            oid = oid_base + obs_index * 3 + (idx % 2)
            obs_lines.append("%d\t%s\t%.17g\t%.17g\t%.17g" %
                             (oid, sid, tick_from_freq(freq, clock), raw.real, raw.imag))
            truth_pred_lines.append("%d\t%.17g\t%.17g" % (oid, la, ph))
            obs_index += 1
    write_lines(os.path.join(data_dir, "observations.tsv"), obs_lines)

    query_freqs = [math.exp(-2.25 + 0.41 * j) for j in range(16)]
    query_lines = ["# query_id station_id instrument_tick"]
    truth_query_lines = []
    for j, freq in enumerate(query_freqs):
        sid = stations[(j * 2 + idx) % len(stations)][0]
        qid = 7000 + idx * 100 + j * 5
        query_lines.append("%d\t%s\t%.17g" % (qid, sid, tick_from_freq(freq, clock)))
        z = forward_model(freq, model, constants)
        la, ph = log_polar(z, constants["pi"])
        truth_query_lines.append("%d\t%.17g\t%.17g" % (qid, la, ph))
    write_lines(os.path.join(data_dir, "queries.tsv"), query_lines)

    write_lines(os.path.join(truth, "truth_clock.tsv"), ["\t".join("%.17g" % x for x in clock)])
    write_lines(os.path.join(truth, "truth_gains.tsv"),
                ["%s\t%.17g\t%.17g\t%.17g\t%.17g" % (sid, gain.real, gain.imag, pole, zero)
                 for sid, _angle, gain, pole, zero in stations])
    write_lines(os.path.join(truth, "truth_model.tsv"), [
        "rho0 %.17g" % model["rho0"], "rho1 %.17g" % model["rho1"], "rho2 %.17g" % model["rho2"],
        "h0 %.17g" % model["h0"], "h1 %.17g" % model["h1"],
        "m %.17g" % model["m"], "tau %.17g" % model["tau"], "c %.17g" % model["c"],
    ])
    write_lines(os.path.join(truth, "truth_predictions.tsv"), truth_pred_lines)
    write_lines(os.path.join(truth, "truth_queries.tsv"), truth_query_lines)
    return label, bundle, truth


def generated_bundles():
    global _GENERATED_ROOT
    if _GENERATED_ROOT is None:
        _GENERATED_ROOT = tempfile.mkdtemp(prefix="colecole_generated_")
        atexit.register(lambda: shutil.rmtree(_GENERATED_ROOT, ignore_errors=True))
        _GENERATED_BUNDLES[:] = [
            generate_bundle(_GENERATED_ROOT, idx, spec)
            for idx, spec in enumerate(generated_case_specs())
        ]
    return list(_GENERATED_BUNDLES)


_GENERATED_BUNDLES = []


def evaluate(bundle, truth_dir):
    sandbox = run_program(bundle)
    try:
        inputs = manifest_paths(bundle)
        state = os.path.join(sandbox, "state")
        required = [os.path.join(state, "clock.tsv"), os.path.join(state, "gains.tsv"), os.path.join(state, "model.tsv"),
                    os.path.join(state, "predictions.tsv"), os.path.join(sandbox, "answers.tsv")]
        if not all(os.path.exists(p) for p in required):
            return {"ok": False, "reason": "missing mandatory state file"}

        truth_model, truth_gains, truth_clock, truth_pred, truth_query = parse_truth(truth_dir)
        clock = {r[0]: float(r[1]) for r in rows(required[0])}
        clock_keys = tuple("a%d" % i for i in range(len(truth_clock)))
        if any(k not in clock or not math.isfinite(clock[k]) for k in clock_keys):
            return {"ok": False, "reason": "incomplete clock state"}
        if max(abs(clock[k] - truth_clock[i]) for i, k in enumerate(clock_keys)) > 2e-8:
            return {"ok": False, "reason": "clock calibration mismatch"}
        got_clock = tuple(clock[k] for k in clock_keys)
        beacon_logs = [math.log(float(r[2])) for r in rows(inputs["clock"])]
        beacon_bounds = (min(beacon_logs), max(beacon_logs))

        gain_rows = rows(required[1])
        got_gains = {
            r[0]: (float(r[1]), float(r[2]), float(r[3]), float(r[4]), float(r[5]))
            for r in gain_rows
        }
        if set(got_gains) != set(truth_gains):
            return {"ok": False, "reason": "station gain coverage mismatch"}
        max_gain = 0.0
        for k, (re, im, pole, zero) in truth_gains.items():
            max_gain = max(max_gain, abs(got_gains[k][0] - re), abs(got_gains[k][1] - im),
                           abs(got_gains[k][2] - pole), abs(got_gains[k][3] - zero))
            if got_gains[k][4] > 1e-6:
                return {"ok": False, "reason": "%s calibration RMS too large" % k}
        if max_gain > GAIN_TOL:
            return {"ok": False, "reason": "gain error %.3g" % max_gain}

        model_rows = rows(required[2])
        model = {r[0]: float(r[1]) for r in model_rows if len(r) >= 2}
        keys = ("rho0", "rho1", "rho2", "h0", "h1", "m", "tau", "c", "objective")
        if any(k not in model or not math.isfinite(model[k]) for k in keys):
            return {"ok": False, "reason": "incomplete model state"}
        if model["objective"] > 0.02:
            return {"ok": False, "reason": "robust objective %.5f too large" % model["objective"]}

        station_angle = {}
        for r in rows(inputs["stations"]):
            station_angle[r[0]] = math.radians(float(r[1]))
        observed = {}
        for r in rows(inputs["observations"]):
            oid, sid = int(r[0]), r[1]
            raw = complex(float(r[3]), float(r[4]))
            gr, gi, pole, zero, _rms = got_gains[sid]
            freq = frequency_from_tick(float(r[2]), got_clock, beacon_bounds)
            z = raw * (1.0 + 1j * 2.0 * math.pi * freq * pole) / (
                complex(gr, gi) * (1.0 + 1j * 2.0 * math.pi * freq * zero)
            ) / cmath.exp(1j * station_angle[sid])
            observed[oid] = (math.log10(abs(z)), math.degrees(cmath.phase(z)))

        expected_obs_ids = [int(r[0]) for r in rows(inputs["observations"])]
        pred_rows = rows(required[3])
        if len(pred_rows) != len(truth_pred) or [int(r[0]) for r in pred_rows] != expected_obs_ids:
            return {"ok": False, "reason": "prediction row count mismatch"}
        max_log = 0.0
        max_phase = 0.0
        for r in pred_rows:
            oid = int(r[0])
            if oid not in truth_pred:
                return {"ok": False, "reason": "unknown observation id"}
            la, ph = float(r[1]), float(r[2])
            if not (math.isfinite(la) and math.isfinite(ph)):
                return {"ok": False, "reason": "non-finite forward state"}
            max_log = max(max_log, abs(la - truth_pred[oid][0]))
            max_phase = max(max_phase, abs(cdiff(ph, truth_pred[oid][1])))
            ola, oph = observed[oid]
            if abs((la - ola) - float(r[3])) > 2e-8 or abs(cdiff(ph, oph) - float(r[4])) > 2e-7:
                return {"ok": False, "reason": "inconsistent residual state"}
        if max_log > LOG_TOL or max_phase > PHASE_TOL:
            return {"ok": False, "reason": "forward state errors %.5f %.4f" % (max_log, max_phase)}

        expected_query_ids = [int(r[0]) for r in rows(inputs["queries"])]
        ans_rows = rows(required[4])
        if [int(r[0]) for r in ans_rows] != expected_query_ids:
            return {"ok": False, "reason": "query order or coverage mismatch"}
        got = {int(r[0]): (float(r[1]), float(r[2])) for r in ans_rows}
        if set(got) != set(truth_query):
            return {"ok": False, "reason": "query coverage mismatch"}
        qlog = qphase = 0.0
        for q, (la, ph) in truth_query.items():
            qlog = max(qlog, abs(got[q][0] - la))
            qphase = max(qphase, abs(cdiff(got[q][1], ph)))
        ok = qlog <= LOG_TOL and qphase <= PHASE_TOL
        return {"ok": ok, "max_log": qlog, "max_phase": qphase,
                "max_pred_log": max_log, "max_pred_phase": max_phase,
                "reason": "" if ok else "query tolerance exceeded"}
    finally:
        shutil.rmtree(sandbox, ignore_errors=True)


def bundles():
    fixed = [
        ("public", os.path.join(APP_DIR, "input_data", "public"), os.path.join(HERE, "reference", "public")),
        ("bundle_a", os.path.join(HERE, "holdout", "bundle_a"), os.path.join(HERE, "reference", "bundle_a")),
        ("bundle_b", os.path.join(HERE, "holdout", "bundle_b"), os.path.join(HERE, "reference", "bundle_b")),
    ]
    return fixed + generated_bundles()
