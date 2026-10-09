"""Advisory benchmark comparisons; timing changes never set a failing exit code."""
import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import random
import statistics

ENV_KEYS = ('python', 'implementation', 'platform', 'architecture', 'cpu',
            'cpu_count', 'host_id', 'affinity', 'six')
METHOD_KEYS = ('timer', 'gc', 'allocation', 'execution')


def load(path, side=None):
    raw = Path(path).read_bytes()
    data = json.loads(raw)
    if not isinstance(data, dict):
        raise ValueError('benchmark artifact must be a JSON object')
    if data.get('schema') == 'gem-benchmark-v1':
        if side:
            raise ValueError('side applies only to historical paired artifacts')
        data['artifact_sha256'] = hashlib.sha256(raw).hexdigest()
        return data
    # Historical definitions were not versioned. Only opposite sides of the
    # exact same paired artifact have shared workload/environment provenance.
    if side not in ('before', 'after') or not isinstance(data.get('rounds'), list):
        raise ValueError('expected gem-benchmark-v1 or a paired artifact with --*-side')
    digest = hashlib.sha256(raw).hexdigest()
    results = {}
    for index, block in enumerate(data['rounds']):
        for name, pair in block.items():
            if side not in pair:
                continue
            entry = results.setdefault(name, {'definition_sha256': digest + ':' + name,
                'unit': 'historical named workload call', 'blocks': []})
            if name in data.get('sample_counts', {}):
                entry['output_points'] = data['sample_counts'][name].get(side)
            entry['blocks'].append(dict(pair[side], pairing_id=digest + ':' + str(index)))
    return {'schema': 'gem-benchmark-v1', 'legacy_artifact': digest,
            'environment': data.get('environment', {}),
            'methodology': data.get('methodology', {}), 'results': results}


def samples(entry):
    values, counts = [], []
    for block in entry['blocks']:
        raw = block['seconds_per_operation']
        if not raw or any(isinstance(v, bool) or not isinstance(v, (int, float)) or
                          not math.isfinite(v) or v <= 0 for v in raw):
            raise ValueError('timing samples must be finite and positive')
        if (isinstance(block['trials'], bool) or not isinstance(block['trials'], int) or
                block['trials'] != len(raw) or isinstance(block['operations_per_trial'], bool) or
                not isinstance(block['operations_per_trial'], int) or block['operations_per_trial'] <= 0):
            raise ValueError('invalid trial/operation counts')
        value = statistics.median(raw) * 1e6
        if not math.isfinite(value):
            raise ValueError('timing magnitude cannot be represented in microseconds')
        values.append(value)
        counts.append(len(raw))
    if not values:
        raise ValueError('empty timing blocks')
    return values, counts


def summary(values):
    median = statistics.median(values)
    return {'median_us': median, 'mad_us': statistics.median(abs(v-median) for v in values),
            'min_us': min(values), 'max_us': max(values), 'block_medians_us': values}


def interval(before, after, paired=False):
    rng = random.Random(305)
    ratios = []
    for _ in range(2000):
        if paired:
            indices = rng.choices(range(len(before)), k=len(before))
            b = [before[i] for i in indices]
            a = [after[i] for i in indices]
        else:
            b = rng.choices(before, k=len(before))
            a = rng.choices(after, k=len(after))
        ratios.append(statistics.median(a) / statistics.median(b))
    ratios.sort()
    return [ratios[49], ratios[1949]]


def compare(baselines, currents, relative_margin=.15, absolute_margin_us=.05):
    if (not baselines or not currents or not math.isfinite(relative_margin) or relative_margin < 0 or
            not math.isfinite(absolute_margin_us) or absolute_margin_us < 0):
        raise ValueError('require measurements and finite nonnegative margins')
    for collection in (baselines, currents):
        digests = [d.get('artifact_sha256', d.get('legacy_artifact')) for d in collection]
        digests = [d for d in digests if d is not None]
        if len(digests) != len(set(digests)):
            raise ValueError('duplicate measurement artifact cannot count as repeated evidence')
    all_data = baselines + currents
    legacy = all_data[0].get('legacy_artifact')
    shared_legacy = legacy is not None and all(d.get('legacy_artifact') == legacy for d in all_data)
    mismatches = []
    first = all_data[0]
    if not shared_legacy:
        for key in ENV_KEYS:
            value = first.get('environment', {}).get(key)
            if value is None or value == 'unavailable' or any(d.get('environment', {}).get(key) != value for d in all_data):
                mismatches.append('environment.' + key)
        for key in METHOD_KEYS:
            value = first.get('methodology', {}).get(key)
            if value is None or any(d.get('methodology', {}).get(key) != value for d in all_data):
                mismatches.append('methodology.' + key)
    output = {'schema': 'gem-comparison-v1', 'advisory_only': True,
        'environment_mismatches': mismatches,
        'environments': [d.get('environment', {}) for d in all_data],
        'baseline_captures': len(baselines), 'current_captures': len(currents),
        'methodologies': [d.get('methodology', {}) for d in all_data],
        'sources': [{'source': d.get('source'),
                     'artifact_sha256': d.get('artifact_sha256', d.get('legacy_artifact'))}
                    for d in all_data],
        'provenance': 'same historical paired artifact' if shared_legacy else 'explicit host and environment metadata',
        'policy': {'relative_review_margin': relative_margin, 'absolute_review_margin_us': absolute_margin_us,
                   'minimum_blocks_per_side': 3, 'minimum_trials_per_block': 5,
                   'bootstrap': '2000 seeded median-ratio resamples; conditional descriptive interval'},
        'results': {}}
    names = sorted(set().union(*(d['results'].keys() for d in all_data)))
    for name in names:
        bentries = [d['results'].get(name) for d in baselines]
        aentries = [d['results'].get(name) for d in currents]
        row = output['results'][name] = {}
        if any(e is None for e in bentries + aentries):
            row['status'] = 'missing_baseline' if any(e is None for e in bentries) else 'missing_current'
            continue
        entries = bentries + aentries
        identities = [(e.get('definition_sha256'), e.get('definition_version'),
                       e.get('metadata'), e.get('unit'), e.get('output_points')) for e in entries]
        if identities[0][0] is None or any(identity != identities[0] for identity in identities):
            row['status'] = 'workload_changed'
            continue
        bv, bc, av, ac = [], [], [], []
        for collection, values, counts in ((bentries, bv, bc), (aentries, av, ac)):
            for entry in collection:
                v, c = samples(entry)
                values.extend(v)
                counts.extend(c)
        row.update(baseline=summary(bv), current=summary(av), unit=entries[0]['unit'])
        ratio = row['current']['median_us'] / row['baseline']['median_us']
        delta = row['current']['median_us'] - row['baseline']['median_us']
        row.update(slowdown_ratio=ratio, speedup_ratio=1/ratio, delta_us=delta)
        if mismatches:
            row['status'] = 'environment_not_comparable'
            continue
        if min(len(bv), len(av)) < 3 or min(bc + ac) < 5:
            row['status'] = 'insufficient_evidence'
            continue
        bids = [b.get('pairing_id') for e in bentries for b in e['blocks']]
        aids = [b.get('pairing_id') for e in aentries for b in e['blocks']]
        paired = shared_legacy and bids == aids and len(set(bids)) == len(bids)
        bounds = interval(bv, av, paired)
        within = []
        for collection in (bentries, aentries):
            deviations = []
            for entry in collection:
                for block in entry['blocks']:
                    raw = block['seconds_per_operation']
                    med = statistics.median(raw)
                    deviations.append(statistics.median(abs(v-med) for v in raw)/med)
            within.append(statistics.median(deviations))
        row['within_block_relative_mad'] = {'baseline': within[0], 'current': within[1]}
        noise = 2 * (max(within[0], row['baseline']['mad_us']/row['baseline']['median_us']) +
                     max(within[1], row['current']['mad_us']/row['current']['median_us']))
        margin = max(relative_margin, noise)
        row.update(ratio_interval_95=bounds, effective_review_margin=margin,
                   resampling='paired blocks' if paired else 'independent blocks')
        regression_votes = sum(v > row['baseline']['median_us']*(1+margin) for v in av)/len(av)
        improvement_votes = sum(v < row['baseline']['median_us']/(1+margin) for v in av)/len(av)
        if bounds[0] > 1+margin and delta > absolute_margin_us and regression_votes >= .8:
            row['status'] = 'possible_regression'
        elif bounds[1] < 1/(1+margin) and -delta > absolute_margin_us and improvement_votes >= .8:
            row['status'] = 'possible_improvement'
        else:
            row['status'] = 'inconclusive'
    output['counts'] = dict(Counter(row['status'] for row in output['results'].values()))
    return output


def markdown(report):
    lines = ['Advisory comparison; flags require review, not automatic acceptance or rejection.', '',
             'Environment differences: ' + (', '.join(report['environment_mismatches']) or 'none recorded'), '']
    for index, env in enumerate(report['environments']):
        side = 'baseline' if index < report['baseline_captures'] else 'current'
        lines.append('%s capture %d: %s; %s; CPU %s; host %s' % (side, index + 1,
            env.get('python', 'unknown'), env.get('platform', 'unknown'),
            env.get('cpu', 'unknown'), env.get('host_id', 'unverified historical artifact')))
    lines.extend(['', '| Workload | Baseline µs | Current µs | Current/baseline | Status |',
                  '|---|---:|---:|---:|---|'])
    for name, row in report['results'].items():
        if 'baseline' in row:
            lines.append('| %s | %.4f | %.4f | %.3f | %s |' % (name,
                row['baseline']['median_us'], row['current']['median_us'], row['slowdown_ratio'], row['status']))
        else:
            lines.append('| %s | — | — | — | %s |' % (name, row['status']))
    return '\n'.join(lines) + '\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=Path, action='append', required=True)
    parser.add_argument('--current', type=Path, action='append', required=True)
    parser.add_argument('--baseline-side', choices=('before', 'after'))
    parser.add_argument('--current-side', choices=('before', 'after'))
    parser.add_argument('--output', type=Path)
    parser.add_argument('--markdown', type=Path)
    parser.add_argument('--relative-margin', type=float, default=.15)
    parser.add_argument('--absolute-margin-us', type=float, default=.05)
    args = parser.parse_args()
    try:
        report = compare([load(p, args.baseline_side) for p in args.baseline],
                         [load(p, args.current_side) for p in args.current],
                         args.relative_margin, args.absolute_margin_us)
    except (ValueError, KeyError, TypeError, OSError) as error:
        parser.error(str(error))
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    rendered = markdown(report)
    if args.markdown:
        args.markdown.parent.mkdir(parents=True, exist_ok=True)
        args.markdown.write_text(rendered)
    print(rendered)


if __name__ == '__main__':
    main()
