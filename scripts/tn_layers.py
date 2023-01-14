import numpy as np
import scanpy as sc
import anndata


def add_layer(output,sample,tn_dict):
    if len(tn_dict)>1:               # 2 layers 2 transposases
        bc_tnh=tn_dict['tnh']
        bc_tn5=tn_dict['tn5']
        adatas_5 = [sc.read(f'{output}/{sample}/{sample}_BC_{bc}.h5ad') for bc in bc_tn5]
        adatas_h = [sc.read(f'{output}/{sample}/{sample}_BC_{bc}.h5ad') for bc in bc_tnh]
        cells = adatas_5[0].obs_names
        for x in range(len(adatas_5)):
            adatas_5[x] = adatas_5[x][cells]
            adatas_h[x] = adatas_h[x][cells]
        A5 = np.sum([adatas_5[x].X for x in range(len(adatas_5))]) 
        AH = np.sum([adatas_h[x].X for x in range(len(adatas_h))]) 
        adata = adatas_5[0].copy()
        adata.X = A5
        adata.layers['tn5'] = A5
        adata.layers['tnH'] = AH
        adata.var['commonness_tn5'] = np.sum(adata.layers['tn5'] >0, axis=0).A1
        adata.var['commonness_tnH'] = np.sum(adata.layers['tnH'] >0, axis=0).A1
        adata.var['commonness'] = adata.var['commonness_tn5'] + adata.var['commonness_tnH']    
        adata.obs['sum_peaks_tn5'] = np.sum(adata.layers['tn5'] >0, axis=1).A1
        adata.obs['coverage_tn5'] = np.sum(adata.layers['tn5'], axis=1).A1
        adata.obs['sum_peaks_tnH'] = np.sum(adata.layers['tnH'] >0, axis=1).A1
        adata.obs['coverage_tnH'] = np.sum(adata.layers['tnH'], axis=1).A1
        adata.obs['sum_peaks'] = adata.obs['sum_peaks_tn5'] + adata.obs['sum_peaks_tnH']
        adata.obs['coverage'] = adata.obs['coverage_tn5'] + adata.obs['coverage_tnH']


        adata.write(f'{output}/{sample}/{sample}.h5ad')
    elif len(tn_dict)==1:             # 1 layer only 1 transposase
        tn=str(tn_dict.keys())[12:15] #name of the transposase used
        bcs=tn_dict[tn]
        adatas = [sc.read(f'{output}/{sample}/{sample}_BC_{bc}.h5ad') for bc in bcs]
        cells = adatas[0].obs_names
        for x in range(4):
            adatas[x] = adatas[x][cells]
        A5 = np.sum([adatas[x].X for x in range(len(adatas))]) 
        adata = adatas[0].copy()
        adata.X = A5
        adata.var['commonness'] = np.sum(adata.X >0, axis=0).A1
        adata.obs['sum_peaks'] = np.sum(adata.X >0, axis=1).A1
        adata.obs['coverage'] = np.sum(adata.X, axis=1).A1
        adata.write(f'{output}/{sample}/{sample}.h5ad')
