#!/usr/bin/env python3
"""Compare independent patient groups using existing Griffin composite coverage TSVs."""
import argparse,csv,hashlib,itertools,json,math
from pathlib import Path
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import f_oneway,t
import statsmodels.formula.api as smf
from statsmodels.stats.anova import anova_lm

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--manifest',required=True,help='TSV: sample, group, path; optional TF columns');ap.add_argument('--output',required=True);ap.add_argument('--subset',default='ALL',help='Named directory under sample directory, e.g. ALL or FILTERED');ap.add_argument('--central-bp',type=int,default=30);ap.add_argument('--plot-bp',type=int,default=1000);ap.add_argument('--permutations',type=int,default=9999);ap.add_argument('--seed',type=int,default=42);ap.add_argument('--covariate',action='append',default=[],help='Optional numeric TSV column; repeatable');a=ap.parse_args()
 mf=Path(a.manifest).resolve();out=Path(a.output).resolve()
 if out.exists() and any(out.iterdir()):raise SystemExit('Output must be new or empty')
 m=pd.read_csv(mf,sep='\t',dtype={'sample':str,'group':str,'path':str});req={'sample','group','path'}
 if not req.issubset(m.columns) or m[list(req)].isna().any().any():raise SystemExit('Manifest requires nonempty sample, group, path')
 if m['sample'].duplicated().any():raise SystemExit('Duplicate sample IDs')
 groups=list(dict.fromkeys(m.group));
 if len(groups) not in (2,3):raise SystemExit('Exactly two or three groups required')
 if min(m.group.value_counts())<2:raise SystemExit('At least two patients per group required')
 if a.central_bp<=0 or a.plot_bp<a.central_bp or a.permutations<1:raise SystemExit('Invalid window or permutations')
 curves=[];positions=None;audit=[]
 for _,r in m.iterrows():
  path=Path(r['path']).expanduser();path=path if path.is_absolute() else mf.parent/path
  if path.is_dir():
   root=path/a.subset if (path/a.subset).is_dir() else path
   cand=sorted(root.rglob('*.GC_corrected.coverage.tsv'))
   if len(cand)!=1:raise SystemExit(f'{r["sample"]}: {len(cand)} candidate TSVs in {root}; provide exact file path or subset')
   path=cand[0]
  d=pd.read_csv(path,sep='\t'); numeric=[]
  for col in d.columns:
   try: numeric.append((int(col),col))
   except (ValueError,TypeError):pass
  numeric.sort();x=np.array([z[0] for z in numeric]); cols=[z[1] for z in numeric]
  if not len(x) or len(set(x))!=len(x):raise SystemExit(f'{path}: missing/duplicate numeric positions')
  if len(d)==1:row=d.iloc[0]
  elif 'site_name' in d.columns and (d.site_name=='COMPOSITE').sum()==1:row=d.loc[d.site_name=='COMPOSITE'].iloc[0]
  else:raise SystemExit(f'{path}: cannot identify exactly one composite row')
  y=pd.to_numeric(row[cols],errors='coerce').to_numpy(float)
  if not np.isfinite(y).all():raise SystemExit(f'{path}: missing or nonfinite positions')
  if x.min()>-a.plot_bp or x.max()<a.plot_bp-5:raise SystemExit(f'{path}: plot window incomplete; examine grid')
  if positions is None:positions=x
  elif not np.array_equal(x,positions):raise SystemExit('Position grids do not match')
  curves.append(y);audit.append({'sample':r['sample'],'group':r['group'],'coverage_path':str(path.resolve()),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
 y=np.stack(curves);x=positions;central=np.abs(x)<=a.central_bp
 if not central.any() or x[central].min()!=-a.central_bp or x[central].max()!=a.central_bp:raise SystemExit('Central boundaries are not represented on grid')
 m['central_coverage']=y[:,central].mean(1);v=m.central_coverage.to_numpy(); idx=[np.flatnonzero(m.group.eq(g)) for g in groups]
 # Permutation ANOVA and full-curve Euclidean PERMANOVA: same label shuffles.
 obs=float(f_oneway(*(v[i] for i in idx)).statistic);mask=np.abs(x)<=a.plot_bp;Y=y[:,mask];D=((Y[:,None,:]-Y[None,:,:])**2).sum(axis=2);n=len(m);k=len(groups);total=D.sum()/(2*n)
 def pseudo(indices):
  within=sum(D[np.ix_(i,i)].sum()/(2*len(i)) for i in indices); between=total-within
  return (between/(k-1))/(within/(n-k)) if within>0 else float('nan')
 fcurve=pseudo(idx);r2=(total-sum(D[np.ix_(i,i)].sum()/(2*len(i)) for i in idx))/total if total else float('nan');rng=np.random.default_rng(a.seed);ex1=ex2=0
 for _ in range(a.permutations):
  perm=rng.permutation(n);pi=[perm[i] for i in idx];ex1+=f_oneway(*(v[i] for i in pi)).statistic>=obs-1e-12;ex2+=pseudo(pi)>=fcurve-1e-12
 rows=[{'test':'Central coverage permutation ANOVA','N':n,'statistic':obs,'P':(ex1+1)/(a.permutations+1),'note':'Global group test; random label permutations'}, {'test':'Full curve Euclidean PERMANOVA','N':n,'statistic':fcurve,'P':(ex2+1)/(a.permutations+1),'R2':r2,'note':f'Positions |x| ≤ {a.plot_bp} bp; random label permutations'}]
 for col in a.covariate:
  if col not in m:raise SystemExit(f'Covariate missing: {col}')
  q=m[['central_coverage','group',col]].copy();q[col]=pd.to_numeric(q[col],errors='coerce');q=q.dropna(); counts=q.group.value_counts()
  if len(counts)!=k or min(counts)<2 or q[col].nunique()<2:rows.append({'test':f'ANCOVA adjusted: {col}','N':len(q),'P':np.nan,'note':'Not estimable: missing groups, too few samples, or constant covariate'});continue
  q=q.rename(columns={col:'covariate'});reduced=smf.ols('central_coverage ~ covariate',data=q).fit();full=smf.ols('central_coverage ~ covariate + C(group)',data=q).fit()
  if full.df_resid<=0 or np.linalg.matrix_rank(full.model.exog)<full.model.exog.shape[1]:rows.append({'test':f'ANCOVA adjusted: {col}','N':len(q),'P':np.nan,'note':'Rank deficient or no residual degrees of freedom'});continue
  test=anova_lm(reduced,full).iloc[1];ranges=q.groupby('group').covariate.agg(['min','max']);overlap=max(ranges['min'])<=min(ranges['max']);rows.append({'test':f'ANCOVA adjusted: {col}','N':len(q),'statistic':test['F'],'P':test['Pr(>F)'],'note':f'Omnibus F; common covariate range across groups: {overlap}; check model assumptions'})
 out.mkdir(parents=True,exist_ok=True);pd.DataFrame(rows).to_csv(out/'tests.tsv',sep='\t',index=False);m.to_csv(out/'patients.tsv',sep='\t',index=False);pd.DataFrame(audit).to_csv(out/'input_audit.tsv',sep='\t',index=False);pd.DataFrame(y,index=m['sample'],columns=x).to_csv(out/'coverage_curves.tsv',sep='\t');m.groupby('group',sort=False).central_coverage.agg(['count','mean','median','std']).to_csv(out/'group_summary.tsv',sep='\t')
 fig,axes=plt.subplots(1,2,figsize=(15,6),gridspec_kw={'width_ratios':[1.5,1]},layout='constrained');colors=['#0072B2','#009E73','#D55E00']
 for j,(g,ii) in enumerate(zip(groups,idx)):
  yy=y[ii][:,mask];avg=yy.mean(0);ci=t.ppf(.95,len(ii)-1)*yy.std(0,ddof=1)/math.sqrt(len(ii));axes[0].plot(x[mask],avg,label=f'{g} (n={len(ii)})',color=colors[j]);axes[0].fill_between(x[mask],avg-ci,avg+ci,color=colors[j],alpha=.15);axes[1].boxplot(v[ii],positions=[j],widths=.4,showfliers=False);axes[1].scatter(j+np.linspace(-.1,.1,len(ii)),v[ii],color=colors[j],s=28,zorder=3)
 axes[0].axvline(0,color='gray',ls=':',lw=.8);axes[0].set(xlabel='Distance from region midpoint (bp)',ylabel='Input composite coverage');axes[0].legend();axes[1].set(xticks=range(k),xticklabels=groups,ylabel=f'Central coverage (±{a.central_bp} bp)');notes='\n'.join(f"{r['test']}: P={r['P']:.4g}" if np.isfinite(r.get('P',np.nan)) else f"{r['test']}: not estimable" for r in rows);axes[1].text(1.04,.98,notes,transform=axes[1].transAxes,va='top',fontsize=9,bbox=dict(facecolor='white',edgecolor='lightgray',boxstyle='round'));fig.suptitle('Griffin composite coverage by group');fig.text(.02,-.02,'Shading: pointwise 90% t confidence intervals; input profiles used as supplied. Independent samples assumed.',fontsize=9);fig.savefig(out/'griffin_group_comparison.png',dpi=220,bbox_inches='tight');fig.savefig(out/'griffin_group_comparison.pdf',bbox_inches='tight');plt.close(fig)
 (out/'analysis_metadata.json').write_text(json.dumps({'subset':a.subset,'central_bp_inclusive':a.central_bp,'plot_bp_inclusive':a.plot_bp,'permutations':a.permutations,'seed':a.seed,'groups':groups,'covariates':a.covariate,'normalization':'input supplied as-is','confidence_interval':'pointwise 90% t','multiple_testing':'not adjusted','permdisp':'not performed','unit':'independent patient/sample'},indent=2)+'\n');print('Saved:',out);print(pd.DataFrame(rows).to_string(index=False))
if __name__=='__main__':main()
