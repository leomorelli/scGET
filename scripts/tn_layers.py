import scanpy as sc
import anndata


def add_layer(output,sample,tn_dict):
    if len(tn_dict)>1:               # 2 layers 2 transposases
        bc_tnh=tn_dict['tnh']
        bc_tn5=tn_dict['tn5']
        adatas_5 = [sc.read(f'{output}/{sample}_BC_{bc}.h5ad') for bc in bc_tn5]
        adatas_h = [sc.read(f'{output}/{sample}_BC_{bc}.h5ad') for bc in bc_tnh]
        adata_5 = adatas_5[0].concatenate(adatas[1:])
        adata_h = adatas_h[0].concatenate(adatas[1:])
        adata_5.layers['tnh']=adata_h.X
        adata_5.layers['tn5']=adata_5.X
        adata_5.write(f'{output}/adata_{sample}.h5ad')
    elif len(tn_dict)==1:             # 1 layer only 1 transposase
        tn=str(tn_dict.keys())[12:15] #name of the transposase used
        bcs=tn_dict[tn]
        adatas = [sc.read(f'{output}/{sample}_BC_{bc}.h5ad') for bc in bcs]
        adata = adatas[0].concatenate(adatas[1:])
        adata.layers[str(tn)]=adata.X
        adata.write(f'{output}/adata_{sample}.h5ad')

