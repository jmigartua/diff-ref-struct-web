#!/usr/bin/env python3
"""Synthetic data files for the formats of which no real example was found (L062), written to manuals/gsas2/synthetic/.

They hold one made-up pattern (a peak on a flat background, 10–30° in steps of 0.02°) in each layout, as GSAS-II's
importers read it (GSASII/imports/G2pwd_*.py, commit fe944cd): PANalytical XRDML, Bruker BRML (zip of XML), Bruker
RAW versions 2 and 4 (binary, offsets from G2pwd_BrukerRAW.py), Rigaku RAS, RASX (zip) and TXT, a CSV table, a Fit2D
.chi and a .qye. They test a reader only against GSAS-II's reading of the same bytes, not against an instrument's file.

Usage (from the repository root): python3 core/tests/data/make_synthetic_formats.py
"""
import io
import math
import os
import struct
import zipfile

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
OUT = os.path.join(ROOT, 'manuals', 'gsas2', 'synthetic')
X = [10.0 + 0.02 * i for i in range(1001)]
Y = [float(round(100 + 1000 * math.exp(-((x - 20.0) / 0.1) ** 2) + 3 * math.sin(7 * x))) for x in X]
FIXED = (2026, 10, 10, 12, 0, 0)  # zip entries' date, so that the files are written the same every time


def zip_bytes(entries):
    b = io.BytesIO()
    with zipfile.ZipFile(b, 'w', zipfile.ZIP_DEFLATED) as z:
        for name, data in entries:
            z.writestr(zipfile.ZipInfo(name, FIXED), data, compress_type=zipfile.ZIP_DEFLATED)
    return b.getvalue()


def xrdml():
    counts = ' '.join(str(int(y)) for y in Y)
    return f'''<?xml version="1.0" encoding="UTF-8"?>
<xrdMeasurements xmlns="http://www.xrdml.com/XRDMeasurement/1.5" status="Completed">
  <comment><entry>synthetic pattern for diff-ref-struct's reader tests</entry></comment>
  <sample type="To be analyzed"><id>synthetic-1</id><name>peak</name></sample>
  <xrdMeasurement measurementType="Scan" status="Completed">
    <usedWavelength intended="K-Alpha 1"><kAlpha1 unit="Angstrom">1.5405980</kAlpha1><kAlpha2 unit="Angstrom">1.5444260</kAlpha2><kBeta unit="Angstrom">1.3922500</kBeta><ratioKAlpha2KAlpha1>0.5000</ratioKAlpha2KAlpha1></usedWavelength>
    <incidentBeamPath><radius unit="mm">240.00</radius><xRayTube id="1" name="Cu LFF"><tension unit="kV">45</tension><current unit="mA">40</current></xRayTube></incidentBeamPath>
    <scan appendNumber="0" mode="Continuous" scanAxis="Gonio" status="Completed">
      <header><startTimeStamp>2026-10-10T12:00:00+02:00</startTimeStamp><endTimeStamp>2026-10-10T12:30:00+02:00</endTimeStamp></header>
      <dataPoints>
        <positions axis="2Theta" unit="deg" listType="StartEnd"><startPosition>10.0000</startPosition><endPosition>30.0000</endPosition></positions>
        <positions axis="Omega" unit="deg" listType="StartEnd"><startPosition>5.0000</startPosition><endPosition>15.0000</endPosition></positions>
        <commonCountingTime unit="seconds">1.00</commonCountingTime>
        <counts unit="counts">{counts}</counts>
      </dataPoints>
    </scan>
  </xrdMeasurement>
</xrdMeasurements>
'''


def brml():
    datums = '\n'.join(f'<Datum>1,1,{x:.4f},{y / 2:.1f},{y:.1f}</Datum>' for x, y in zip(X, Y))
    raw = f'''<?xml version="1.0" encoding="utf-8"?>
<RawData xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <TimeStampStarted>2026-10-10T12:00:00</TimeStampStarted>
  <TimeStampFinished>2026-10-10T12:30:00</TimeStampFinished>
  <Identifier><SampleName>synthetic-1</SampleName></Identifier>
  <DataRoutes>
    <DataRoute><ScanInformation><MeasurementPoints>3</MeasurementPoints></ScanInformation></DataRoute>
    <DataRoute>
      <ScanInformation><MeasurementPoints>{len(X)}</MeasurementPoints><TimePerStepEffective>1</TimePerStepEffective></ScanInformation>
{datums}
    </DataRoute>
  </DataRoutes>
</RawData>
'''
    return zip_bytes([('Experiment0/RawData0.xml', raw)])


def ras():
    rows = '\n'.join(f'{x:.4f} {y:.0f} 1.0000' for x, y in zip(X, Y))
    return ('*RAS_DATA_START\n*RAS_HEADER_START\n*FILE_SAMPLE "synthetic-1"\n*HW_XG_WAVE_LENGTH_ALPHA1 "1.540593"\n'
            '*HW_XG_WAVE_LENGTH_ALPHA2 "1.544414"\n*MEAS_SCAN_STEP "0.0200"\n*RAS_HEADER_END\n*RAS_INT_START\n'
            f'{rows}\n*RAS_INT_END\n*RAS_DATA_END\n')


def rasx():
    prof = '﻿' + ''.join(f'{x:.4f}\t{y:.0f}\t1\r\n' for x, y in zip(X, Y))
    cond = ('<?xml version="1.0" encoding="utf-8"?><MeasurementConditions><SampleName>synthetic-1</SampleName>'
            '<WavelengthKalpha1>1.540593</WavelengthKalpha1><WavelengthKalpha2>1.544414</WavelengthKalpha2></MeasurementConditions>')
    return zip_bytes([('Data0/Profile0.txt', prof.encode('utf-8')), ('Data0/MesurementConditions0.xml', cond)])


def rigaku_txt():
    out = ['synthetic pattern', 'exported as text', '']
    for i in range(0, 1000, 4):
        out.append(f'{X[i]:.4f} ' + ' '.join(f'{y:.0f}' for y in Y[i:i + 4]))
    return '\n'.join(out) + '\n'


def csv():
    return '# 2Theta, Intensity, ESD\n' + ''.join(f'{x:.4f},{y:.1f},{math.sqrt(y):.4f}\n' for x, y in zip(X, Y))


def chi():
    return ('synthetic pattern\n2-Theta Angle (Degrees)\nIntensity\n' + f'{len(X):12d}\n' +
            ''.join(f'  {x:.6E}  {y:.6E}\n' for x, y in zip(X, Y)))


def qye():
    lam = 0.4136
    rows = ''.join(f'{4 * math.pi * math.sin(math.radians(x / 2)) / lam:.6f} {y:.1f} {math.sqrt(y):.4f}\n' for x, y in zip(X, Y))
    return f'# synthetic pattern in Q\n# wavelength {lam}\n' + rows


def raw_v2():
    head = bytearray(256)
    head[0:4] = b'RAW2'
    struct.pack_into('<i', head, 4, 1)
    head[168:188] = b'10/10/26 12:00:00'.ljust(20, b'\0')
    head[188:190] = b'Cu'
    struct.pack_into('<fff', head, 190, 1.540598, 1.544426, 0.5)
    struct.pack_into('<f', head, 206, 1.39222)
    block = bytearray(48)
    struct.pack_into('<HH', block, 0, 48, len(X))
    struct.pack_into('<ff', block, 12, 0.02, 10.0)
    return bytes(head) + bytes(block) + struct.pack(f'<{len(Y)}f', *Y)


def raw_v4():
    b = bytearray(61)
    b[0:7] = b'RAW4.00'
    b[12:24] = b'10/10/2026'.ljust(12, b'\0')
    b[24:34] = b'12:00:00'.ljust(10, b'\0')
    # A key/value segment (type 10), the X-ray source (type 30), then the range (160).
    kv = struct.pack('<II', 10, 8 + 4 + 24 + 16) + b'\0' * 4 + b'Operator'.ljust(24, b'\0') + b'synthetic'.ljust(16, b'\0')
    src = struct.pack('<II', 30, 120) + b'\0' * 64 + struct.pack('<ddddd', 1.5418, 1.540598, 1.544426, 1.39222, 0.5) + b'\0' * 4 + b'Cu'.ljust(4, b'\0')
    rng = bytearray(struct.pack('<I', 160))
    # The range header as GSAS-II reads it from just after the 160 mark: scan type at +28, start, step, points and
    # step time at +68, generator at +96, wavelength at +108, datum and header sizes at +132; header segments at +156.
    body = bytearray(156)
    body[28:52] = b'Locked Coupled'.ljust(24, b'\0')
    struct.pack_into('<ddIf', body, 68, 10.0, 0.02, len(X), 1000.0)
    struct.pack_into('<ff', body, 96, 40.0, 40.0)
    struct.pack_into('<d', body, 108, 1.540598)
    struct.pack_into('<II', body, 132, 4, 72)
    seg = struct.pack('<II', 50, 72) + b'\0' * 4 + b'2Theta'.ljust(24, b'\0') + b'\0' * 20 + struct.pack('<d', 10.0) + b'\0' * 8
    return bytes(b) + kv + src + bytes(rng) + bytes(body) + seg + struct.pack(f'<{len(Y)}f', *Y)


def main():
    os.makedirs(OUT, exist_ok=True)
    files = {
        'scan.xrdml': xrdml(), 'scan.brml': brml(), 'scan.ras': ras(), 'scan.rasx': rasx(), 'scan.txt': rigaku_txt(),
        'scan.csv': csv(), 'scan.chi': chi(), 'scan.qye': qye(), 'scan_v2.raw': raw_v2(), 'scan_v4.raw': raw_v4(),
    }
    for name, data in files.items():
        with open(os.path.join(OUT, name), 'wb') as f:
            f.write(data.encode('utf-8') if isinstance(data, str) else data)


if __name__ == '__main__':
    main()
