import argparse
from argparse import RawTextHelpFormatter
import sys
import pandas as pd

parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="""
=================================================
 Filter PSI output based on min/max sample PSI
=================================================
""")
parser.add_argument('--i', action='store',
                    required=True,
                    help='Input PSI file')
parser.add_argument('--o', action='store',
                    required=True,
                    help='Output filtered PSI file')
args = parser.parse_args()

if args.i == args.o:
    sys.exit('Output filename can not be the same as input')

df_main = pd.read_csv(args.i, sep='\t')

slope = 0.5
maxPSIboundary = 25
dropList = []
count = 0
for ix in df_main.index:
    count = count + 1
    psi = df_main.loc[ix].tolist()[df_main.columns.get_loc('ExonBoundary')+1:]
    psi = [float(x) for x in psi]
    psi_max = max(psi)
    if psi_max < 0:#maxPSIboundary:
        dropList.append(ix)
    else:
        psi_min = min([x for x in psi if x >= 0])
        dmm = psi_max - psi_min
        y = slope*psi_max + maxPSIboundary - 0.25*maxPSIboundary
        if dmm < y:
            dropList.append(ix)
    if count % 1000 == 0:
        print(f'{count} exons processed from {args.i}')
print(f'{count} total exons processed')
df_main.drop(dropList, inplace=True)

df_main.to_csv(args.o, sep='\t', index=False)
