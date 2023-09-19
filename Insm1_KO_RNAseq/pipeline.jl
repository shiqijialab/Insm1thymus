include("../deps/dataProcessKit/dataProcessKit.jl")
using .dataProcessKit

## E18 samples
#{{ Mapping
KK(path="STAR_Mapping")do
    smp=readtb("../sample_info_E18.tsv")
    smp["cmd"]=map(smp["fileID"], smp["sample"]) do fn, lb
        "STAR --runThreadN 8 --readFilesCommand zcat --sjdbGTFfile /path/to/Mus_musculus.GRCm38.92.addChr.gtf --genomeDir /path/to/STAR_index/mm10.92.chr.150nt --outFilterMismatchNmax 3 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --readFilesIn ../mouse_RNAseq/fq_gz/Mouse_X150-01-$(fn)_good_1.fq.gz ../mouse_RNAseq/fq_gz/Mouse_X150-01-$(fn)_good_2.fq.gz --outFileNamePrefix $(lb)- --outSAMattributes AS NH HI nM MD"
    end
    qsub_sh(smp["cmd"], name="STAR", pernode=4, qsubparam="-lnodes=1:ppn=36", tryone=true)
end

#Sort bam file by readID
KK(path="_Bam_sortedByName")do
    smp=jload("../sample_info")
    smp["cmd"]=map(smp["sample"]) do lb
        "samtools sort -n -T tmp4sort.$lb -o $lb.sortedByName.bam ../STAR_Mapping/$(lb)-Aligned.sortedByCoord.out.bam"
    end
    qsub_sh(smp["cmd"], name="sortbam", pernode=20, qsubparam="-lnodes=1:ppn=20")
end
#}}
#{{ HTSeq counting
KK(path="HTSeq")do
    smp=readtb("../sample_info_E18.tsv")
    C=map(smp["sample"]) do lb
        """
htseq-count -f bam -r name -s no ../_Bam_sortedByName/$lb.sortedByName.bam /path/to/Mus_musculus.GRCm38.92.addChr.gtf > $lb-htseq-nostrand.tsv.tmp
mv $lb-htseq-nostrand.tsv.tmp $lb-htseq-nostrand.tsv
"""
    end
    qsub_sh(C, name="thymus_htseq", pernode=20, qsubparam="-lnodes=1:ppn=20")
end
## Read number
smp=jload("sample_info")
C=map(smp["sample"]) do lb
    T=readtb("HTSeq/$lb-htseq-nostrand.tsv", head=c"geneid, N::readnum")
    rec(T, .!startswith.(T["geneid"], "_"))
end
@assert iselequal(x->x["geneid"], C)
ge=tb(geneid=C[1]["geneid"], readnum=hcat(map(x->x["readnum"], C)...))
jsave("gene_readnum", ge, meta=ds(smpID=smp["sample"]))
#}}
#{{ Gene expression
ge, smpID=jload2("gene_readnum", "smpID")
smp=readtb("../sample_info_E18.tsv")
@assert smpID==smp["sample"]
geinf=jload("/path/to/Mus_musculus.GRCm38.79.gene.jld2")
ge=tbx(ge, geinf, "geneid", auniq=true, buniq=true)
ge["FPKM"]=(1e+9)*ge["readnum"]./sum(ge["readnum"], dims=1)./ge["gene_len"]

# Add TAR info
t=readlines("/path/to/TRA_gene_list/Aire_dependent_TRAs_geneid.txt")
ge["AireDepTRA"]=ismbr(ge["geneid"], t)
t=readlines("/path/to/TRA_gene_list/Aire_independent_TRAs_geneid.txt")
ge["AireIndepTRA"]=ismbr(ge["geneid"],t)
    
jsave("gene_readnum_FPKM", ge, meta=ds(smpID=smp["sample"]))
#}}
#{{ DESeq2 test
KK()do
    smp=readtb("../sample_info_E18.tsv")
    ge, smpID=jload2("gene_readnum_FPKM", "smpID")
    @assert smp["sample"]==smpID
    ge=rec(ge,sumh(ge["readnum"]).>=2*rnum(smp))
    t=DESeq2(ge["readnum"][:, smp["treatment"].=="WT"], ge["readnum"][:, smp["treatment"].=="M"])
    ge["raw_pvalue"]=t["pvalue"]
    ge["padj"]=t["padj"]
    ge["logFC"]=t["log2FoldChange"]
    jsave("E18_gene_DESeq2_test", ge, meta=ds(smpID=smpID))
end
#}}
#{{ Volcano plot
using myDraw
ge, smp=jload2("E18_gene_DESeq2_test.jld2", "smpID")
figure(figsize=(4, 4))
yy=min.(-log10.(ge["raw_pvalue"]), 10)
xx=ge["logFC"]
grp=@withrow ge begin
    if ($padj<0.1) .& (abs.($logFC).>log2(1.5))
        $logFC<0 ? 4 : 5
    elseif ($raw_pvalue<0.05) .& (abs.($logFC).>log2(1.5))
        $logFC<0 ? 2 : 3
    else
        1
    end
end

colorls=[blendcolor("w", "b", 0.667), blendcolor("w", "r"), "b", "r"]
grploop(fastgrp(grp, 1:5), xx, yy, id=true) do gi, x, y
    if gi==1
        kdescatter(x, y, s=2, cmap="binary", normto=(0.3, 1), rasterized=true)
    else
        plot(x, y, ".", ms=3, color=colorls[gi-1], rasterized=true)
    end
    #plot(x, y, ".", ms=3, color=["0.75", colorls...][gi])
end
N=freqas(2:5, grp)
for (x, y, n, c) in zip([0.05, 0.95, 0.05, 0.95], [3, 3, 10, 10], N, colorls)
    txtplotrl(x, y, string(n), color=c, rl=:x)
end

fdr01=@with ge (minimum($raw_pvalue[$padj.>0.1])+maximum($raw_pvalue[$padj.<0.1]))/2
guideline("-", -log10(fdr01), label="FDR=0.1")
guideline("-", -log10(0.05), label="P=0.05")
xlabel(k"$log_2$ fold change")
ylabel(k"$-log_{10}(P)$")
ticklabelreplace(:y, "10"=>"≥10")
l=ismbr(ge["gene_name"], c"Aire, Insm1")
txtplot(xx[l], yy[l], ge["gene_name"][l], repel=true)

pltsave("Volcano_plot_E18_Insm1_ko2")
#}}
#{{ Mini heatmap
using myDraw
ge, smpID=jload2("E18_gene_DESeq2_test.jld2", "smpID")
smp=jload("../Insm1_KO_megSortCell_RNAseq/sample_info.jld2")
@assert smp["sample"]==smpID
rename!(ge, "logFC"=>"lgFC", "padj"=>"FDR", "raw_pvalue"=>"rawP")
ge["isTRA"]=ge["AireDepTRA"] .| ge["AireIndepTRA"]
delete!(ge, "AireDepTRA")
delete!(ge, "AireIndepTRA")

gA, smpA=jload2("Aire_KO_RNAseq_DESeq2.jld2", "smp")
gA=rec(gA, (gA["FDR"].<0.1) .& (gA["lgFC"].<-1))
ge["isAireDep"]=ismbr(ge["geneid"], gA["geneid"])
ge["isAireIndep"]=ge["isTRA"] .& .!ge["isAireDep"]

e=rec(ge, (ge["isTRA"] .& (ge["rawP"].<0.05) .& (ge["lgFC"].<-log2(1.5))) .| ismbr(ge["gene_name"], c"Insm1, Aire"))
ge=rec(ge, sortri((sign.(ge["lgFC"]),ge["gene_name"])))
ge["zFPKM"]=zscore(log2.(ge["FPKM"][:, [1,3,5,2,4,6]].+1), 2)

M=fill(0.0, 5+count(ge["isTRA"]), 6)
offset=0.000000001
M[1,:]=ge["zFPKM"][ge["gene_name"].=="Insm1",:].+offset
M[3,:]=ge["zFPKM"][ge["gene_name"].=="Aire",:].+offset
t=count(ge["isAireDep"])
M[5:5+t-1, :]=ge["zFPKM"][ge["isAireDep"],:].+offset
M[5+t+1:end, :]=ge["zFPKM"][ge["isAireIndep"],:].+offset

M=[M[:,1:3] fill(0.0,size(M,1)) M[:,4:6]]
genels=[c"Insm1,,Aire,";ge["gene_name"][ge["isAireDep"]];"";ge["gene_name"][ge["isAireIndep"]]]

figure(figsize=(3,12))
xlim(1,size(M,2)+1)
ylim(1,size(M,1)+1)
ax=gca()
ax[:invert_yaxis]()
for i=1:size(M,1), j=1:size(M,2)
    v=M[i,j]
    v==0 && continue
    ax[:add_patch](rectangle((j,i),1,1,fc=cmap(clrfun(v)),ec="k", linewidth=0.5))
end
halfbox(left=false,bottom=false)
xticks([1.5:3.5; 5.5:7.5], 
    c"Wt1, Wt2, Wt3, M1, M2, M3")
yticks(1.5:(size(M,1)+0.5), genels, fontsize=7);
pltsave("figures/mini_heatmap_Insm1ko_RNAseq_E18")
#}}

## Adult samples
#{{ Mapping
KK(path="STAR_Mapping")do
    smp=readtb("../sample_info_adult.tsv")
    smp["cmd"]=map(smp["label"]) do lb
        "STAR --runThreadN 16 --readFilesCommand zcat --sjdbGTFfile /path/to/Mus_musculus.GRCm38.92.addChr.gtf --genomeDir /path/to/STAR_index/mm10.92.chr.150nt --outFilterMismatchNmax 3 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --readFilesIn ../fastq/Mouse_AE732-01-$(lb)_good_1.fq.gz ../fastq/Mouse_AE732-01-$(lb)_good_2.fq.gz --outFileNamePrefix $(lb)- --outSAMattributes AS NH HI nM MD"
    end
    qsub_sh(smp["cmd"], name="STAR")
end

#Sort bam file by readID
KK(path="_Bam_sortedByName")do
    smp=readtb("../sample_info.tsv")
    smp["cmd"]=map(smp["label"]) do lb
        "samtools sort -n -T tmp4sort.$lb -o $lb.sortedByName.bam ../STAR_Mapping/$(lb)-Aligned.sortedByCoord.out.bam"
    end
    qsub_sh(smp["cmd"], name="sortbam", pernode=16)
end
#}}
#{{ HTSeq counting
KK(path="HTSeq")do
    smp=readtb("../sample_info_adult.tsv")
    C=map(smp["label"]) do lb
        """
htseq-count -f bam -r name -s no ../_Bam_sortedByName/$lb.sortedByName.bam /path/to/Mus_musculus.GRCm38.79.addChr.gtf > $lb-htseq-nostrand.tsv.tmp
mv $lb-htseq-nostrand.tsv.tmp $lb-htseq-nostrand.tsv
"""
    end
    qsub_sh(C, name="thymus_htseq", pernode=16)
end
## Read number
smp=readtb("sample_info.tsv")
C=map(smp["label"]) do lb
    T=readtb("HTSeq/$lb-htseq-nostrand.tsv", head=c"geneid, N::readnum")
    rec(T, .!startswith.(T["geneid"], "_"))
end
@assert iselequal(x->x["geneid"], C)
ge=tb(geneid=C[1]["geneid"], readnum=hcat(map(x->x["readnum"], C)...))
jsave("gene_readnum", ge, meta=ds(smpID=smp["label"]))
#}}
#{{ Gene expression
ge, smpID=jload2("gene_readnum", "smpID")
smp=readtb("../sample_info_adult.tsv")
@assert smpID==smp["label"]
geinf=jload("/path/to/Mus_musculus.GRCm38.79.gene.jld2")
ge=tbx(ge, geinf, "geneid", auniq=true, buniq=true)
ge["FPKM"]=(1e+9)*ge["readnum"]./sum(ge["readnum"], dims=1)./ge["gene_len"]

# Add TAR info
t=readlines("/path/to/TRA_gene_list/Aire_dependent_TRAs_geneid.txt")
ge["AireDepTRA"]=ismbr(ge["geneid"], t)
t=readlines("/path/to/TRA_gene_list/Aire_independent_TRAs_geneid.txt")
ge["AireIndepTRA"]=ismbr(ge["geneid"],t)
    
jsave("gene_readnum_FPKM", ge, meta=ds(smpID=smp["label"]))
#}}
#{{ DESeq2 test
KK()do
    smp=readtb("../sample_info_adult.tsv")
    ge, sample_name=jload2("adult_Insm1KO_megBatch1_gene_readnum_FPKM", "sample_name")
    @assert smp["sample_name"]==sample_name
    
    ge=rec(ge, meanh(ge["readnum"]).>=2)
    t=DESeq2(ge["readnum"][:, smp["treatment"].=="WT"],
             ge["readnum"][:, smp["treatment"].=="M"])
    ge["rawP"]=t["pvalue"]
    ge["FDR"]=t["padj"]
    ge["lgFC"]=t["log2FoldChange"]
    
    jsave("adult_gene_DESeq2_test", ge, meta=ds(sample_name=sample_name))
end
#}}
#{{ Volcano plot
using myDraw
ge, smp=jload2("adult_gene_DESeq2_test.jld2", "sample_name")
figure(figsize=(4, 4))
yy=-log10.(ge["rawP"])
xx=ge["lgFC"]
grp=@withrow ge begin
    if ($FDR<0.1) .& (abs.($lgFC).>log2(1.5))
        $lgFC<0 ? 4 : 5
    elseif ($rawP<0.05) .& (abs.($lgFC).>log2(1.5))
        $lgFC<0 ? 2 : 3
    else
        1
    end
end
colorls=[blendcolor("w", "b", 0.667), blendcolor("w", "r"), "b", "r"]
grploop(fastgrp(grp, 1:5), xx, yy, id=true) do gi, x, y
    if gi==1
        kdescatter(x, y, s=2, cmap="binary", normto=(0.3, 1), rasterized=true)
    else
        plot(x, y, ".", ms=3, color=colorls[gi-1], rasterized=true)
    end
#    plot(x, y, ".", ms=3, color=["0.75", colorls...][gi])
end
N=freqas(2:5, grp)
for (x, y, n, c) in zip([0.05, 0.95, 0.05, 0.95], [3, 3, 9, 9], N, colorls)
    txtplotrl(x, y, string(n), color=c, rl=:x)
end

fdr01=@with ge (minimum($rawP[$FDR.>0.1])+maximum($rawP[$FDR.<0.1]))/2
guideline("-", -log10(fdr01), label="FDR=0.1")
guideline("-", -log10(0.05), label="P=0.05")
#guideline("|", log2(1.5).*[1, -1])
#guideline("|", "k:")
xlabel(k"$log_2$ fold change")
ylabel(k"$-log_{10}(P)$")

l=ismbr(ge["gene_name"], c"Aire, Insm1")
txtplot(xx[l], yy[l], ge["gene_name"][l], repel=true)

pltsave("Volcano_plot_adult_Insm1_ko2")
#}}
#{{ Mini heatmpa
using myDraw
ge, smpID=jload2("adult_gene_DESeq2_test.jld2", "sample_name")
smp=jload("../Insm1_KO_megSortCell_RNAseq/sample_info.jld2")
@assert smp["sample"]==smpID
ge["isTRA"]=ge["AireDepTRA"] .| ge["AireIndepTRA"]
delete!(ge, "AireDepTRA")
delete!(ge, "AireIndepTRA")

gA, smpA=jload2("Aire_KO_RNAseq_DESeq2.jld2", "smp")
gA=rec(gA, (gA["FDR"].<0.1) .& (gA["lgFC"].<-1))
ge["isAireDep"]=ismbr(ge["geneid"], gA["geneid"])
ge["isAireIndep"]=ge["isTRA"] .& .!ge["isAireDep"]

e=rec(ge, (ge["isTRA"] .& (ge["rawP"].<0.05) .& (ge["lgFC"].<-log2(1.5))) .| ismbr(ge["gene_name"], c"Insm1, Aire"))
ge=rec(ge, sortri((sign.(ge["lgFC"]),ge["gene_name"])))
ge["zFPKM"]=zscore(log2.(ge["FPKM"][:, [1,3,5,2,4,6]].+1), 2)

M=fill(0.0, 5+count(ge["isTRA"]), 6)
offset=0.000000001
M[1,:]=ge["zFPKM"][ge["gene_name"].=="Insm1",:].+offset
M[3,:]=ge["zFPKM"][ge["gene_name"].=="Aire",:].+offset
t=count(ge["isAireDep"])
M[5:5+t-1, :]=ge["zFPKM"][ge["isAireDep"],:].+offset
M[5+t+1:end, :]=ge["zFPKM"][ge["isAireIndep"],:].+offset

M=[M[:,1:3] fill(0.0,size(M,1)) M[:,4:6]]
genels=[c"Insm1,,Aire,";ge["gene_name"][ge["isAireDep"]];"";ge["gene_name"][ge["isAireIndep"]]]

figure(figsize=(3,12))
xlim(1,size(M,2)+1)
ylim(1,size(M,1)+1)
ax=gca()
ax[:invert_yaxis]()
for i=1:size(M,1), j=1:size(M,2)
    v=M[i,j]
    v==0 && continue
    ax[:add_patch](rectangle((j,i),1,1,fc=cmap(clrfun(v)),ec="k", linewidth=0.5))
end
halfbox(left=false,bottom=false)
xticks([1.5:3.5; 5.5:7.5], 
    c"Wt1, Wt2, Wt3, M1, M2, M3")
yticks(1.5:(size(M,1)+0.5), genels, fontsize=7);
pltsave("figures/mini_heatmap_Insm1ko_RNAseq_adult")
#}}

## Analysis public Aire-KO and Fezf2-KO data
# Raw files were downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144877
#{{ Get Aire data
KK()do
    T=readtb(shc"zcat GSE144877_Aire_KO_RNA_seq.csv.gz", quotes=true, autotype=true, dm=',')
    t=c"""
ctl - R_2017_01_26_18_26_12_user_Immunol-30-170126_takaba_stephanie.IonXpress_001 (single) trimmed (GE) - Gene name=>geneid
Feature ID=>gene_name
EDGE test: ctl vs AireKO, tagwise dispersions - FDR p-value correction=>EDGE_FDR
EDGE test: ctl vs AireKO, tagwise dispersions - Fold change=>EDGE_FC
EDGE test: ctl vs AireKO, tagwise dispersions - P-value=>EDGE_rawP
"""c
    ge=tb()
    for x in t
        (x1, x2)=split(x, "=>")
        if x2 in c"gene_name, geneid"
            ge[x2]=T[x1]
        else
            ge[x2]=parse.(Float64, T[x1])
        end
    end

    smpls=c"""AireKO - R_2017_01_26_18_26_12_user_Immunol-30-170126_takaba_stephanie.IonXpress_004 (single) trimmed (GE)
AireKO - R_2017_01_26_18_26_12_user_Immunol-30-170126_takaba_stephanie.IonXpress_005 (single) trimmed (GE)
AireKO - R_2017_01_26_18_26_12_user_Immunol-30-170126_takaba_stephanie.IonXpress_006 (single) trimmed (GE)
ctl - R_2017_01_26_18_26_12_user_Immunol-30-170126_takaba_stephanie.IonXpress_001 (single) trimmed (GE)
ctl - R_2017_01_26_18_26_12_user_Immunol-30-170126_takaba_stephanie.IonXpress_002 (single) trimmed (GE)
ctl - R_2017_01_26_18_26_12_user_Immunol-30-170126_takaba_stephanie.IonXpress_003 (single) trimmed (GE)
"""c
    @assert iselequal(dict2tuple(T, f"$1 - Gene name".(smpls)))
    ge["readnum"]=hcat(int.(dict2tuple(T, f"$1 - Unique gene reads".(smpls)))...)
    ge["RPKM"]=hcat(float.(dict2tuple(T, f"$1 - RPKM".(smpls)))...)

    #Filter by gene type
    geinf=jload("/home/luna.kuleuven.be/u0121733/projects/biodata/annotation/mm10/Mus_musculus.GRCm38.92.gene.jld2")
    ge["gene_typ"]=dtshift(ge["geneid"], geinf["geneid"], geinf["gene_typ"], "", safe=true)
    ge=rec(ge, ismbr(ge["gene_typ"], c"protein_coding, lincRNA"))

    #DESeq2 test
    ge=@rec(ge, meanh($readnum).>=2)
    T=DESeq2(ge["readnum"][:, 4:6], ge["readnum"][:, 1:3])
    d"FDR, rawP, lgFC"ge=d"padj, pvalue, log2FoldChange"T

    jsave("Aire_KO_RNAseq_DESeq2", ge, meta=ds(smp=smpls))
end
#}}
#{{ Get Fezf2 data
KK()do
    T=readtb(shc"zcat GSE144877_Fezf2_KO_RNA_seq.csv.gz", quotes=true, autotype=true, dm=',')
    t=c"""
Fezf2cKO - R_2017_01_26_18_26_12_user_Immunol-30-170126_takaba_stephanie.IonXpress_010 (single) trimmed (GE) - Gene name=>geneid
Feature ID=>gene_name
EDGE test: ctl vs Fezf2cKO, tagwise dispersions - FDR p-value correction=>EDGE_FDR
EDGE test: ctl vs Fezf2cKO, tagwise dispersions - Fold change=>EDGE_FC
EDGE test: ctl vs Fezf2cKO, tagwise dispersions - P-value=>EDGE_rawP
"""c
    ge=tb()
    for x in t
        (x1, x2)=split(x, "=>")
        if x2 in c"gene_name, geneid"
            ge[x2]=T[x1]
        else
            ge[x2]=parse.(Float64, T[x1])
        end
    end

    smpls=c"""
Fezf2cKO - R_2017_01_26_18_26_12_user_Immunol-30-170126_takaba_stephanie.IonXpress_010 (single) trimmed (GE)
Fezf2cKO - R_2017_01_26_18_26_12_user_Immunol-30-170126_takaba_stephanie.IonXpress_011 (single) trimmed (GE)
Fezf2cKO - R_2017_01_26_18_26_12_user_Immunol-30-170126_takaba_stephanie.IonXpress_012 (single) trimmed (GE)
ctl - R_2017_01_26_18_26_12_user_Immunol-30-170126_takaba_stephanie.IonXpress_007 (single) trimmed (GE)
ctl - R_2017_01_26_18_26_12_user_Immunol-30-170126_takaba_stephanie.IonXpress_008 (single) trimmed (GE)
ctl - R_2017_01_26_18_26_12_user_Immunol-30-170126_takaba_stephanie.IonXpress_009 (single) trimmed (GE)
"""c
    @assert iselequal(dict2tuple(T, f"$1 - Gene name".(smpls)))
    ge["readnum"]=hcat(int.(dict2tuple(T, f"$1 - Unique gene reads".(smpls)))...)
    ge["RPKM"]=hcat(float.(dict2tuple(T, f"$1 - RPKM".(smpls)))...)

    #Filter by gene type
    geinf=jload("/home/luna.kuleuven.be/u0121733/projects/biodata/annotation/mm10/Mus_musculus.GRCm38.92.gene.jld2")
    ge["gene_typ"]=dtshift(ge["geneid"], geinf["geneid"], geinf["gene_typ"], "", safe=true)
    ge=rec(ge, ismbr(ge["gene_typ"], c"protein_coding, lincRNA"))

    #DESeq2 test
    ge=@rec(ge, meanh($readnum).>=2)
    T=DESeq2(ge["readnum"][:, 4:6], ge["readnum"][:, 1:3])
    d"FDR, rawP, lgFC"ge=d"padj, pvalue, log2FoldChange"T

    jsave("Fezf2_KO_RNAseq_DESeq2", ge, meta=ds(smp=smpls))
end
#}}
