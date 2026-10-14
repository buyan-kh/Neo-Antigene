import csv
import json
import os
import re

allele_id_map = {
    'MHC-011': 'HLA-A*02:01',
    'MHC-018': 'HLA-A*01:01',
    'MHC-025': 'HLA-A*03:01',
    'MHC-032': 'HLA-A*11:01',
    'MHC-039': 'HLA-A*24:02',
    'MHC-046': 'HLA-B*07:02',
    'MHC-053': 'HLA-B*08:01',
    'MHC-060': 'HLA-B*44:02',
    'MHC-067': 'HLA-C*07:01',
    'MHC-074': 'HLA-C*08:02',
}

def parse_allele(a_str):
    s = a_str.strip().strip('"\'').strip()
    if s in allele_id_map:
        return allele_id_map[s]
    s_upper = s.upper()
    if s_upper.startswith("HLA-"):
        s_upper = s_upper[4:]
    elif s_upper.startswith("HLA"):
        s_upper = s_upper[3:]
    m = re.match(r'^([ABC])\*?(\d{2}):?(\d{2})$', s_upper)
    if m:
        locus, major, minor = m.groups()
        return f"HLA-{locus}*{major}:{minor}"
    raise ValueError(f"Unknown allele format: {a_str}")

def parse_peptide(p_str):
    return p_str.strip().strip('"\'').strip().upper()

def load_training_pairs():
    training_pairs = set()
    
    # panel_a.csv
    with open("predictor_training_set/panel_a.csv") as f:
        reader = csv.DictReader(f)
        for row in reader:
            p = parse_peptide(row['peptide'])
            a = parse_allele(row['mhc_allele'])
            training_pairs.add((p, a))
            
    # panel_b.tsv
    with open("predictor_training_set/panel_b.tsv") as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            p = parse_peptide(row['sequence'])
            a = parse_allele(row['hla'])
            training_pairs.add((p, a))
            
    # panel_c.csv
    with open("predictor_training_set/panel_c.csv") as f:
        reader = csv.DictReader(f)
        for row in reader:
            p = parse_peptide(row['peptide'])
            a = parse_allele(row['allele_id'])
            training_pairs.add((p, a))
            
    return training_pairs

def compute_auroc(y_true, y_score):
    pos_scores = [s for y, s in zip(y_true, y_score) if y == 1]
    neg_scores = [s for y, s in zip(y_true, y_score) if y == 0]
    if not pos_scores or not neg_scores:
        return None
    pairs_won = 0.0
    for p in pos_scores:
        for n in neg_scores:
            if p > n:
                pairs_won += 1.0
            elif p == n:
                pairs_won += 0.5
    return pairs_won / (len(pos_scores) * len(neg_scores))

def main():
    training_pairs = load_training_pairs()
    
    benchmark_data = []
    with open("benchmark.csv") as f_bm, open("predictions.csv") as f_pred:
        r_bm = csv.DictReader(f_bm)
        r_pred = csv.DictReader(f_pred)
        for row_bm, row_pred in zip(r_bm, r_pred):
            pep = parse_peptide(row_bm['peptide'])
            allele = parse_allele(row_bm['mhc_allele'])
            label = row_bm['true_label'].strip()
            score = float(row_pred['pred_score'])
            
            is_leaked = (pep, allele) in training_pairs
            
            benchmark_data.append({
                'peptide': pep,
                'allele': allele,
                'is_binder': 1 if label == 'binder' else 0,
                'score': score,
                'is_leaked': is_leaked
            })
            
    leaked_pairs = [[d['peptide'], d['allele']] for d in benchmark_data if d['is_leaked']]
    clean_examples = [d for d in benchmark_data if not d['is_leaked']]
    
    clean_y = [d['is_binder'] for d in clean_examples]
    clean_scores = [d['score'] for d in clean_examples]
    clean_auroc = compute_auroc(clean_y, clean_scores)
    
    os.makedirs("output", exist_ok=True)
    
    submission = {
        "contaminated": True,
        "leaked_pair_count": len(leaked_pairs),
        "leaked_pairs": leaked_pairs,
        "mechanism": "train_test_overlap",
        "trusted_metric": False,
        "corrected_metric": {
            "name": "auroc",
            "value": clean_auroc,
            "n_clean_pairs": len(clean_examples),
            "uncomputable": False
        }
    }
    
    with open("output/submission.json", "w") as f:
        json.dump(submission, f, indent=2)
        
    with open("output/leaked_pairs.csv", "w", newline="") as f:
        writer = csv.writer(f)
        for pep, allele in leaked_pairs:
            writer.writerow([pep, allele])
            
    print(f"Generated submission with {len(leaked_pairs)} leaked pairs.")
    print(f"Clean pairs: {len(clean_examples)}, AUROC: {clean_auroc}")

if __name__ == "__main__":
    main()
