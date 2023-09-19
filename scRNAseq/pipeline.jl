include("../deps/dataProcessKit/dataProcessKit.jl")
using .dataProcessKit

#{{ Cell type identificaton
using PyCall
sc=importpy("scanpy")

AD=sc.read_10x_mtx("raw_data", prefix="GSM5831744_adult_perinate_gex_", gex_only=false)
sc.external.pp.scrublet(AD)

#Select the cells according to the barcode list from Michelson et al, Cell 2022.
smpbc=readtb(`zcat raw_data/GSM5831744_adult_perinate_metadata.txt.gz`, dm=" ")
AD=pyf"$1[$2, :]"(AD, ismbr(AD.obs.index.to_list(), smpbc["barcode"]))

set!(AD.obs, "sample", dtshift(AD.obs.index.to_list(), smpbc["barcode"], smpbc["hash_ident"], safe=true))

AD.obs.__setitem__("n_counts", AD.X.sum(1))
AD.obs.__setitem__("n_genes", pyf"($1>0).sum(1)"(AD.X))
ismtgene=startswith.(AD.var.index.to_list(), "mt-")
AD.obs.__setitem__("mt_prop", pyf"($1[:, $2]).sum(1)"(AD.X, ismtgene)./AD.obs.n_counts)

markergenes=c"""
Foxn1, Trp63, Aire, Pou2f3, Spib, Sox8, Hnf4g, Hnf4a,
Foxj1, Trp73, Myog, Foxi1, Foxi2, Ptf1a, Foxa2, Foxa3, Foxa1,
Spdef, Grhl1, Gata3
"""
markergenes2=c"""
Psmb11, Ly75, Mki67, Top2a, Ccl21a, Krt5, Aire, H2-Aa, Chat, Trpm5, Ccl9, Tnfaip2, Apoa4, Guca2a, Dynlrb2, Dnah12,
Actc1, Ckm, Cftr, Atp6v1g3, Prss2, Clps, Chgb, Ret, Aqp4, Muc5b, Spink5, Rptn, Krt10, Krtdap, Ivl, Sbsn
"""
AD=pyf"$1[:, $2]"(AD, (h2v(pyf"($1>0).sum(0)"(AD.X)).>=200) .|
                  ismbr(AD.var.index.to_list(), [markergenes; markergenes2]));

sc.pp.normalize_total(AD, target_sum=1e5)
sc.pp.log1p(AD)
sc.pp.scale(AD)

sc.pp.highly_variable_genes(AD, n_top_genes=2000)

sc.pp.pca(AD, n_comps=50, use_highly_variable=true, svd_solver="arpack")
sc.pl.pca_variance_ratio(AD, n_pcs=50)

sc.pp.neighbors(AD, n_pcs=30)
sc.tl.umap(AD)
sc.pl.umap(AD, color="sample")

sc.tl.leiden(AD, resolution=1.9)
sc.pl.umap(AD, color="leiden")

sc.tl.dendrogram(AD, groupby="leiden")
sc.pl.dotplot(AD, markergenes, groupby="leiden", dendrogram=true)

sc.tl.dendrogram(AD, groupby="leiden")
sc.pl.dotplot(AD, markergenes2, groupby="leiden", dendrogram=true)

t=c"""Perinatal cTEC: 18
Adult cTEC: 16
Transit-amplifying: 14, 15
Immature mTEC: 2, 4, 1, 6
Aire-expressing: 12
Tuft_1: 0, 3
Tuft_2: 5, 11
Microfold: 8
Enterocyte/hepatocyte: 17, 9
Ciliated: 20
Muscle: 19
Neuroendocrine: 10
Basal(lung): 7
Basal(skin)/Keratinocyte: 13
"""c
T=mapr(t, mrows=true) do x
    C=split(x, ':')
    celltype=string(strip(C[1]))
    clusterNO=string.(strip.(split(C[2], ',')))
    tbexp((clusterNO=clusterNO,), celltype=celltype)
end
set!(AD.obs, "celltype", dtshift(AD.obs.leiden.to_list(), T["clusterNO"], T["celltype"], safe=true))
sc.pl.umap(AD, color="celltype")

AD.write("clustered_anndata_resolution1d9_cellAssigned.h5ad")
#}}

#{{ Plot figures
using PyCall
sc=importpy("scanpy")

AD=sc.read_10x_mtx("raw_data", prefix="GSM5831744_adult_perinate_gex_", gex_only=false)
sc.external.pp.scrublet(AD)

smpbc=readtb(`zcat raw_data/GSM5831744_adult_perinate_metadata.txt.gz`, dm=" ")
AD=pyf"$1[$2, :]"(AD, ismbr(AD.obs.index.to_list(), smpbc["barcode"]))

set!(AD.obs, "sample", dtshift(AD.obs.index.to_list(), smpbc["barcode"], smpbc["hash_ident"], safe=true))

AD.obs.__setitem__("n_counts", AD.X.sum(1))
AD.obs.__setitem__("n_genes", pyf"($1>0).sum(1)"(AD.X))
ismtgene=startswith.(AD.var.index.to_list(), "mt-")
AD.obs.__setitem__("mt_prop", pyf"($1[:, $2]).sum(1)"(AD.X, ismtgene)./AD.obs.n_counts)

sc.pp.normalize_total(AD, target_sum=1e5)
sc.pp.log1p(AD)
sc.pp.scale(AD)

AD0=sc.read("clustered_anndata_resolution1d9_cellAssigned.h5ad")
foo=(geneid, score_name; ensembl::Bool=true)->begin
    genenm=if ensembl
        dtshift(geneid, cut".1".(AD.var.gene_ids.to_list()), AD.var.index.to_list(), "", safe=true)
    else
        geneid
    end
    genenm=genenm[.!isempty.(genenm)]
    sc.tl.score_genes(AD, genenm, score_name=score_name)
    @assert pyf"$1.to_list()==$2.to_list()"(AD.obs.index, AD0.obs.index)
    set!(AD0.obs, score_name, AD.obs[score_name])
end

#Adult Insm1-ko
ge=jload("../Insm1_KO_RNAseq/adult_gene_DESeq2_test.jld2", showinfo=false)
gl=@with ge $geneid[($rawP.<0.05) .& ($lgFC.<-log2(1.5)) .& ($AireDepTRA .| $AireIndepTRA)]
foo(gl, f"Insm1KO_TRA_dw_adult_$1genes"(length(gl)))

#E18 Insm1-ko
ge=jload("../Insm1_KO_RNAseq/E18_gene_DESeq2_test.jld2", showinfo=false)
gl=@with ge $geneid[($raw_pvalue.<0.05) .& ($logFC.<-log2(1.5)) .& ($AireDepTRA .| $AireIndepTRA)]
foo(gl, f"Insm1KO_TRA_dw_E18_dw_$1genes"(length(gl)))

ADa=pyf"$1[$2, :]"(AD0, startswith.(get(AD0.obs, "sample").to_list(), "adult"))
ADe=pyf"$1[$2, :]"(AD0, startswith.(get(AD0.obs, "sample").to_list(), "neonate"))

orderls=c"""Neuroendocrine
Enterocyte/hepatocyte
Microfold
Ciliated
Basal(skin)/Keratinocyte
Basal(lung)
Muscle
Tuft_1
Tuft_2
Transit-amplifying
Immature mTEC
Aire-expressing
Perinatal cTEC
Adult cTEC
"""c

sc.pl.violin(ADa, keys="Insm1KO_TRA_dw_adult_81genes", groupby="celltype", rotation=90, 
    order=orderls[orderls.!="Perinatal cTEC"])
gcf().set_size_inches(5.5, 3.5)
txtplotrl(0.9, 0.95, "adults")
pltsave("Insm1KO_TRA_dwInAdult_violin_adultCells")

sc.pl.violin(ADe, keys="Insm1KO_TRA_dw_E18_dw_79genes", groupby="celltype", rotation=90, 
    order=orderls[orderls.!="Adult cTEC"])
gcf().set_size_inches(5.5, 3.5)
txtplotrl(0.1, 0.95, "neonates")
pltsave("Insm1KO_TRA_dwInE18_violin_neonateCells")

orderls_a=orderls[orderls.!="Perinatal cTEC"]
sc.pl.violin(ADa, keys=c"Insm1, Aire", groupby="celltype", ncols=2, rotation=90,
    order=orderls[orderls.!="Perinatal cTEC"])
gcf().set_size_inches(11, 3.5)
pltsave("Insm1_and_Aire_exp_violin_adultCells_authorCellAnn")

sc.pl.violin(ADe, keys=c"Insm1, Aire", groupby="celltype", ncols=2, rotation=90,
    order=orderls_e)
gcf().set_size_inches(11, 3.5)
pltsave("Insm1_and_Aire_exp_violin_neonateCells_authorCellAnn")

markergenes=c"""
Insm1,
Foxn1, Trp63, Aire, Pou2f3, Spib, Sox8, Hnf4g, Hnf4a,
Foxj1, Trp73, Myog, Foxi1, Foxi2, Ptf1a, Foxa2, Foxa3, Foxa1,
Spdef, Grhl1, Gata3
"""
sc.tl.dendrogram(AD, groupby="celltype")
sc.pl.dotplot(AD, markergenes, groupby="celltype", categories_order=orderls, figsize=(7.6, 4))
pltsave("cellmarker_exp_per_celltype_dotPlot_authorCellAnn")

#Load UMAP coordinates provided by Michelson et al, Cell 2022.
authorinf=readtb("raw_data/authors_celltype_umap_coordinates.tsv", autotype=true)
T=dtshift(AD.obs.index.to_list(), authorinf["barcode"], authorinf, safe=true)
set!(AD.obsm, "X_umap", hcat(T["umap1"], T["umap2"]))
sc.pl.umap(AD, color="celltype")
pltsave("UMAP_celltype")

sc.pl.umap(AD, color=c"Insm1, Aire")
pltsave("UMAP_Insm1_and_Aire")
#}}
