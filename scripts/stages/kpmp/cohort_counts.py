import pandas as pd
t=pd.read_csv('kpmp_open_file_manifest.tsv',sep='\t',dtype={'participant_id':str})
t['WSI']=t.wsi_total>0; t['PDS']=t.tiv_descriptor_scores
t['scsn']=(t.snRNA_expr_matrix+t.scRNA_expr_matrix+t.multiome_expr_matrix)>0
t['Visium']=t.visium_expr_matrix>0; t['SEG']=t.seg_masks_files>0
t['core4']=(t.wsi_HE+t.wsi_frozenHE>0)&(t.wsi_PAS>0)&(t.wsi_TRI>0)&(t.wsi_SIL>0)
enr,adj,dm,ht=t['Enrollment Category'],t['Primary Adjudicated Category'].fillna(''),t['Diabetes History'].fillna(''),t['Hypertension History'].fillna('')
groups={
 'DKD_adjudicated': adj=='Diabetic Kidney Disease',
 'CKD_diabetesYes': (enr=='CKD')&(dm=='Yes'),
 'HKD_adjudicated': adj=='Hypertensive Kidney Disease',
 'HKD_adjudicated_DMno': (adj=='Hypertensive Kidney Disease')&(dm=='No'),
 'CKD_DMno_HTNyes': (enr=='CKD')&(dm=='No')&(ht=='Yes'),
 'DM-R': enr=='DM-R','HealthyRef': enr=='Healthy Reference',
 'CKD_all': enr=='CKD','AKI_all': enr=='AKI','ALL': enr.notna()}
rows=[]
for g,msk in groups.items():
    x=t[msk]; w=x.WSI&x.PDS
    rows.append({'group':g,'n':len(x),'WSI':int(x.WSI.sum()),'WSI_HE_PAS_TRI_SIL':int(x.core4.sum()),'TIV_scores':int(x.PDS.sum()),
      'SegMasks':int(x.SEG.sum()),'WSI+TIV':int(w.sum()),'WSI+TIV+sc/sn':int((w&x.scsn).sum()),'WSI+TIV+Visium':int((w&x.Visium).sum()),
      'WSI+TIV+(sc/sn|Vis)':int((w&(x.scsn|x.Visium)).sum()),'WSI+sc/sn':int((x.WSI&x.scsn).sum()),'WSI+Visium':int((x.WSI&x.Visium).sum()),
      'WSI+SEG+TIV':int((w&x.SEG).sum()),'any_sc/sn':int(x.scsn.sum()),'any_Visium':int(x.Visium.sum())})
r=pd.DataFrame(rows); pd.set_option('display.width',300); print(r.to_string(index=False))
r.to_csv('kpmp_cohort_counts.tsv',sep='\t',index=False)
print('\nmedian WSI/participant:',t.loc[t.WSI,'wsi_total'].median(),'; total open WSI TB:',round(t.wsi_total_GB.sum()/1000,2),
      '; DKD WSI TB:',round(t.loc[groups['DKD_adjudicated'],'wsi_total_GB'].sum()/1000,2))
print('in manifest but not clinical csv:',(~t.in_clinical_csv).sum())
print(pd.crosstab(t.loc[t.PDS,'Enrollment Category'],t.loc[t.PDS,'Primary Adjudicated Category'].fillna('NA')))
