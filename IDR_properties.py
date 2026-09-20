import pandas as pd
import numpy as np
from tqdm import tqdm

from localcider.sequenceParameters import SequenceParameters

INPUT_FILE = "disordered_regions.xlsx"
OUTPUT_FILE = "IDR_properties.xlsx"

DISORDER_PROMOTING = set("AEGKPQRS")
ORDER_PROMOTING = set("CFILVWY")

STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")

df = pd.read_excel(INPUT_FILE)

FCR = []
NCPR = []
SCD = []
KAPPA = []
HYDROPATHY = []

POSITIVE = []
NEGATIVE = []

NET_CHARGE = []
CHARGE_DENSITY = []

DISORDER_PERCENT = []
ORDER_PERCENT = []

ERRORS = []

# for every IDR

for seq in tqdm(df["Disorder Region"]):

    # Clean sequence
    seq = str(seq).upper().strip()
    seq = "".join([aa for aa in seq if aa in STANDARD_AA])

    length = len(seq)

    if length == 0:

        FCR.append(np.nan)
        NCPR.append(np.nan)
        SCD.append(np.nan)
        KAPPA.append(np.nan)
        HYDROPATHY.append(np.nan)
       
        POSITIVE.append(np.nan)
        NEGATIVE.append(np.nan)

        NET_CHARGE.append(np.nan)
        CHARGE_DENSITY.append(np.nan)

        DISORDER_PERCENT.append(np.nan)
        ORDER_PERCENT.append(np.nan)

        ERRORS.append("Empty sequence")

        continue

 
    pos = seq.count("K") + seq.count("R")
    neg = seq.count("D") + seq.count("E")

    net = pos - neg

    disorder_count = sum(aa in DISORDER_PROMOTING for aa in seq)
    order_count = sum(aa in ORDER_PROMOTING for aa in seq)

    POSITIVE.append(pos)
    NEGATIVE.append(neg)

    NET_CHARGE.append(net)
    CHARGE_DENSITY.append(round(net / length, 4))

    DISORDER_PERCENT.append(round(disorder_count / length * 100, 2))
    ORDER_PERCENT.append(round(order_count / length * 100, 2))

   
    # localCIDER based
    

    try:

        sp = SequenceParameters(seq)

        FCR.append(round(sp.get_FCR(), 4))
        NCPR.append(round(sp.get_NCPR(), 4))
        SCD.append(round(sp.get_SCD(), 4))
        KAPPA.append(round(sp.get_kappa(), 4))
        HYDROPATHY.append(round(sp.get_mean_hydropathy(), 4))

        ERRORS.append("")

    except Exception as e:

        FCR.append(np.nan)
        NCPR.append(np.nan)
        SCD.append(np.nan)
        KAPPA.append(np.nan)
        HYDROPATHY.append(np.nan)

        ERRORS.append(str(e))


df["FCR"] = FCR
df["NCPR"] = NCPR
df["SCD"] = SCD
df["Kappa"] = KAPPA
df["Mean Hydropathy"] = HYDROPATHY


df["Positive Residues"] = POSITIVE
df["Negative Residues"] = NEGATIVE

df["Net Charge"] = NET_CHARGE
df["Charge Density"] = CHARGE_DENSITY

df["Disorder-promoting Residues (%)"] = DISORDER_PERCENT
df["Order-promoting Residues (%)"] = ORDER_PERCENT

df["Calculation Error"] = ERRORS


df.to_excel(OUTPUT_FILE, index=False)

print("\nDone!")
print(f"Saved to: {OUTPUT_FILE}")
print(f"Total IDRs processed: {len(df)}")