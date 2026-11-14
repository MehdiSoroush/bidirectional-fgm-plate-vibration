# Abaqus 2021 Python 2.7: abaqus python export_u3_abaqus2021.py Job.odb 2 output.csv
from __future__ import print_function
import sys, csv
from odbAccess import openOdb
if len(sys.argv) != 4:
    raise SystemExit('Usage: abaqus python export_u3_abaqus2021.py JOB.odb MODE_NUMBER output.csv')
odb_path, mode_text, out_path = sys.argv[1:]
mode = int(mode_text)
odb = openOdb(path=odb_path, readOnly=True)
try:
    step = list(odb.steps.values())[-1]
    frames = step.frames
    candidates = [f for f in frames if getattr(f, 'mode', None) == mode]
    if len(candidates) != 1:
        raise RuntimeError('Could not uniquely find mode %d; available: %s' % (mode, [(i, getattr(f, 'mode', None), f.description) for i, f in enumerate(frames)]))
    frame = candidates[0]
    if 'U' not in frame.fieldOutputs:
        raise RuntimeError('U field missing; request displacement output for frequency step.')
    coords = {}
    for instance_name, instance in odb.rootAssembly.instances.items():
        for node in instance.nodes:
            coords[(instance_name, node.label)] = node.coordinates
    rows = []
    for v in frame.fieldOutputs['U'].values:
        key = (v.instance.name, v.nodeLabel)
        if key not in coords:
            continue
        xyz = coords[key]
        rows.append((v.instance.name, v.nodeLabel, xyz[0], xyz[1], xyz[2], v.data[2]))
    rows.sort(key=lambda r: (r[0], r[1]))
    with open(out_path, 'wb') as f:
        writer = csv.writer(f)
        writer.writerow(['instance','node','x','y','z','u3'])
        writer.writerows(rows)
    print('Exported %d nodal U3 values from mode %d to %s' % (len(rows), mode, out_path))
finally:
    odb.close()
