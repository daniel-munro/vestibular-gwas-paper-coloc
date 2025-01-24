import argparse
import pandas as pd
from pathlib import Path

def ld_blocks(chrom: list, pos: list, ld_file: Path):
    """Return LD block for each SNP"""
    ld = pd.read_csv(ld_file, sep='\s+')
    ld['loc'] = [f'Loc{n}' for n in range(1, len(ld) + 1)]
    # Split ld into a df per chromosome:
    ld = {chr: ld[ld['chr'] == chr] for chr in ld['chr'].unique()}
    # Expand first and last block to include all SNPs:
    for c in ld.keys():
        ld[c]['start'].iloc[0] = 0
        ld[c]['stop'].iloc[-1] = 1e12
    # Print progress:
    i = 0    
    for c, p in zip(chrom, pos):
        i += 1
        if i % 100000 == 0:
            print(f'Processed {i} of {len(chrom)}')
        ld_blocks = ld[c][(ld[c]['start'] <= p) & (ld[c]['stop'] > p)]['loc']
        if len(ld_blocks) == 0:
            print(c, p)
            yield "NA"
        block = ld_blocks.values[0]
        yield block

parser = argparse.ArgumentParser(description='Convert GWAS summary stats to z-scores in torus input format')
parser.add_argument('--gwas', help='GWAS summary stats file', type=Path, required=True)
parser.add_argument('--format', help='Format of GWAS summary stats file', choices=['mvp', 'gtex'], required=True)
parser.add_argument('--ld', help='LD block BED file', type=Path, required=True)
parser.add_argument('--out', help='Output file', type=Path, required=True)
args = parser.parse_args()

if args.format == 'mvp':
    df = pd.read_csv(args.gwas, sep='\s+', usecols=['SNP', 'CHR', 'BP', 'Beta', 'SE'])
    # df = df.iloc[:100, ]
    df['CHR'] = 'chr' + df['CHR'].astype(str)
    df['LOC'] = list(ld_blocks(df['CHR'], df['BP'], args.ld))
    df['Z'] = (df['Beta'] / df['SE']).round(6)
    df = df[['SNP', 'LOC', 'Z']]
    df.to_csv(args.out, sep='\t', index=False, header=False)
elif args.format == 'gtex':
    df = pd.read_csv(args.gwas, sep='\t', usecols=['variant_id', 'chromosome', 'position', 'zscore'])
    df = df[~df['variant_id'].isna()]
    # df = df.iloc[:100, ]
    df['loc'] = list(ld_blocks(df['chromosome'], df['position'], args.ld))
    df['zscore'] = df['zscore'].round(6)
    df = df[['variant_id', 'loc', 'zscore']]
    df.to_csv(args.out, sep='\t', index=False, header=False)
