import scanpy as sc
import anndata


def add_layer(output,sample,tn_dict):
    if len(tn_dict)>1:               # 2 layers 2 transposases
        bc_tnh=tn_dict['tnh']
        bc_tn5=tn_dict['tn5']
        adatas_5 = [sc.read(f'{output}/{sample}/{sample}_BC_{bc}.h5ad') for bc in bc_tn5]
        adatas_h = [sc.read(f'{output}/{sample}/{sample}_BC_{bc}.h5ad') for bc in bc_tnh]
        adata0_5 = adatas_5.pop(0)
        adata0_h = adatas_h.pop(0)
        cells = adata0_5.obs_names
        for ad5 in adatas_5:
            adata0_5.X = adata0_5.X + ad5[cells].X
            adata0_5.obs['n_regions'] += ad5[cells].obs['n_regions']
            adata0_5.var['n_cells'] += ad5[cells].var['n_cells']
        for adh in adatas_h:
            adata0_h.X = adata0_h.X + adh[cells].X
            adata0_h.obs['n_regions'] += adh[cells].obs['n_regions']
            adata0_h.var['n_cells'] += adh[cells].var['n_cells']
        adata0_5.layers['tnh']=adata0_h.X
        adata0_5.layers['tn5']=adata0_5.X
        adata0_5.write(f'{output}/{sample}/{sample}.h5ad')
    elif len(tn_dict)==1:             # 1 layer only 1 transposase
        tn=str(tn_dict.keys())[12:15] #name of the transposase used
        bcs=tn_dict[tn]
        adatas = [sc.read(f'{output}/{sample}/{sample}_BC_{bc}.h5ad') for bc in bcs]
        adata0 = adatas.pop(0)
        cells = adata0.obs_names
        for ad in adatas:
            adata0.X = adata0.X + ad[cells].X
            adata0.obs['n_regions'] += ad[cells].obs['n_regions']
            adata0.var['n_cells'] += ad[cells].var['n_cells']
        adata0.write(f'{output}/{sample}/{sample}.h5ad')
