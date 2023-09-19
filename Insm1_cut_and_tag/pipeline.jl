include("../deps/dataProcessKit/dataProcessKit.jl")
using .dataProcessKit

#{{ Adopter trimming
KK(path="_fastq_trimmed") do
    smp=readtb("../sample_info.tsv")
    C=map(smp["label"], smp["fastq1"], smp["fastq2"]) do lb, fq1, fq2
        "fastp --in1 $fq1 --in2 $fq2 --out1 $(lb)_R1.trimmed.fq.gz --out2 $(lb)_R2.trimmed.fq.gz"
    end
    bgrunsh(C, name="trim")
end
#}}
#{{ Bowtie2 Mapping
KK(path="_Bowtie2_Mapping") do
    smp=readtb("../sample_info.tsv")
    C=map(smp["label"]) do lb
        "bowtie2 -p 4 --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -X 700 -x /path/to/Bowtie_index/mm10/mm10_ -1 ../_fastq_trimmed/$(lb)_R1.trimmed.fq.gz -2 ../_fastq_trimmed/$(lb)_R2.trimmed.fq.gz -S $lb.sam"
    end
    bgrunsh(C, name=f"map_$1".(smp["label"]))
end
## Handle sam file
KK(cdto="_Bowtie2_Mapping") do
    smp=readtb("../sample_info.tsv")
    C=map(smp["label"]) do lb
        "samtools view -b $lb.sam | samtools sort -T _tmp_$(lb) > $lb.sorted.bam"
    end
    bgrunsh(C, name="map")
end
## Convert to bigwig format
KK(cdto="_Bowtie2_Mapping") do
    smp=readtb("../sample_info.tsv")
    C=map(smp["label"]) do lb
        "bamCoverage -b $lb.sorted.bam -o $lb.bw"
    end
    bgrunsh(C, name="bigwig")
end
#}}
#{{ Peak calling
KK(path="Peak_calling")do
    smp=readtb("../sample_info.tsv")
    C, _=grpfun(smp["individual"], smp) do T
        pos=T["label"][findonly(.!T["isNeg"])]
        neg=T["label"][findonly(T["isNeg"])]
        "macs2 callpeak -t ../_Bowtie2_Mapping/$pos.sorted.bam -c ../_Bowtie2_Mapping/$neg.sorted.bam -n autodup_$(pos)_$(neg) -g mm --keep-dup auto"
    end
    bgrunsh(C, name="pkcll")
end
## Parse all peakcalls
KK(cdto="Peak_calling") do
    fns=matchfile("*_peaks.narrowPeak", onlyone=false)
    lbs=cut"_1:end-1".(fns)
    for lb in lbs
        T=readtb(sfc"grep -v '^#' $(1)_peaks.xls | sed '/^$/d'"(lb), autotype=true)
        T["chrno"]=chr2no(T["chr"])
        T=rec(T, T["chrno"].>0)
        T["pos"]=[T["start"] T["end"]]
        delete!.((T,), c"start, end")
        jsave("$(lb)_peaks", T)
    end
end
#}}
#{{ Merge E18 Insm1 peak
Insm1pks=readlines("E18_Insm1_cutAndTag_peak_file_list.txt")
T=vcatr_with("smpNO"=>1:3, jload.(Insm1pks)...)
T=rec(T, T["-log10(qvalue)"].>5)
T["pkNO"]=1:rnum(T)

foo=(A, B)->begin
    Ta=rec(T, T["smpNO"].==A)
    Tb=rec(T, T["smpNO"].==B)
    ai, bi=genomemap(d"chr, pos"Ta, d"chr, pos"Tb, touch=true)
    [Ta["pkNO"][ai] Tb["pkNO"][bi]]
end
M=vcat(foo(1, 2), foo(2, 3), foo(1, 3))
G=linkgrp(M)
pk, _ = grpfun(G, M) do x
    I=unival(vec(x))
    Is=I[argmax(T["-log10(qvalue)"][I])]
    smpNO=join(string.(unique(T["smpNO"][I]))) 
    
    tb(chr=T["chr"][Is],
       chrno=T["chrno"][Is],
       pos=[minimum(T["pos"][I, 1]) maximum(T["pos"][I, 2])],
       abs_summit=T["abs_summit"][Is],
       max_score=T["-log10(qvalue)"][Is],
       smpNO=smpNO)
end
ti, pi=genomemap(d"chrno, pos"T, d"chrno, pos"pk, touch=true)
T["max_score"]=T["-log10(qvalue)"]
T["smpNO"]=string.(T["smpNO"])
pk=vcatr(pk, rec_i(fd(T, keys(pk)), ti))

jsave("merged_Insm1_peak_E18", pk)

## Prepare sequences for motif scan
pk=jload("merged_Insm1_peak_E18")
seq=getseq(d"chrno, pos"pk, path="/path/to/genome/mm10/mm10_jld")
writefasta("Insm1_loose_peak_seq.fa", seq, f"$1:$2-$3".(pk["chr"], pk["pos"][:, 1], pk["pos"][:, 2]))
#}}
#{{ Merge adult Insm1 peak
Insm1pks=readlines("Adult_Insm1_cutAndTag_peak_file_list.txt")
T=vcatr_with("smpNO"=>string.(1:2), jload.(Insm1pks)...)
T=rec(T, T["-log10(qvalue)"].>5)
T["pkNO"]=1:rnum(T)

foo=(A, B)->begin
    Ta=rec(T, T["smpNO"].==A)
    Tb=rec(T, T["smpNO"].==B)
    ai, bi=genomemap(d"chr, pos"Ta, d"chr, pos"Tb, touch=true)
    [Ta["pkNO"][ai] Tb["pkNO"][bi]]
end
M=foo("1", "2")
G=linkgrp(M)
pk, _ = grpfun(G, M) do x
    I=unival(vec(x))
    Is=I[argmax(T["-log10(qvalue)"][I])]
    smpNO=if all(T["smpNO"][I].=="1")
        "1"
    elseif all(T["smpNO"][I].=="2")
        "2"
    else
        "12"
    end     
    tb(chr=T["chr"][Is],
       chrno=T["chrno"][Is],
       pos=[minimum(T["pos"][I, 1]) maximum(T["pos"][I, 2])],
       abs_summit=T["abs_summit"][Is],
       max_score=T["-log10(qvalue)"][Is],
       smpNO=smpNO)
end
ti, pi=genomemap(d"chrno, pos"T, d"chrno, pos"pk, touch=true)
T["max_score"]=T["-log10(qvalue)"]
pk=vcatr(pk, rec_i(fd(T, keys(pk)), ti))
jsave("merged_Insm1_peak_Adult", pk)
#}}
#{{ Insm1 peak binding genes
geinf=jload("/home/luna.kuleuven.be/u0121733/projects/biodata/annotation/mm10/Mus_musculus.GRCm38.92.gene.jld2")
geinf["scan_pos"]=rowfun(geinf["gene_pos"], geinf["strand"]) do p, s
    if s=='+'
        [p[1]-2000 p[1]+500]
    else
        @assert s=='-'
        [p[2]-500 p[2]+2000]
    end
end
pk=jload("merged_Insm1_peak_E18.jld2")
# pk=jload("merged_Insm1_peak_Adult.jld2")
gi, pi=genomemap(d"chrno, scan_pos"geinf, d"chrno, pos"pk, touch=true)
geinf["isInsm1"]=falsesbut(rnum(geinf), gi)

jsave("Insm1_promoter_binding_gene_in_E18", geinf)
# jsave("Insm1_promoter_binding_gene_in_Adult", geinf)
#}}
#{{ Insm1 peak in different gene regions
KK()do
    exn=jload("/path/to/mm10/Mus_musculus.GRCm38.92.jld2")
    cdg=rec(exn, ismbr(exn["geneid"], exn["geneid"][exn["cds_pos"][:, 1].>0]))

    utr5=rowpile()
    utr3=rowpile()
    cdg=rec(cdg, sortri(d"transid, exn_pos"cdg))
    grploop(cdg["transid"], cdg["exn_pos"], cdg["cds_pos"], cdg["strand"], cdg["chrno"]) do ep, cp, strand, chrno
        utrL=rowmx()
        for i=1:size(ep, 1)
            if cp[i, 1]<=0
                addrow!(utrL, ep[i, :])
            else
                ep[i, 1] < cp[i, 1]-1 && addrow!(utrL, [ep[i, 1] cp[i, 1]-1])
                break
            end
        end
        utrL=something(value(utrL), zerorow(Int, 2))

        utrR=rowmx()
        for i=size(ep, 1):-1:1
            if cp[i, 1]<=0
                addrow!(utrR, ep[i, :])
            else
                cp[i, 2]+1<ep[i, 2] && addrow!(utrR, [cp[i, 2]+1 ep[i, 2]])
                break
            end
        end
        utrR=something(value(utrR), zerorow(Int, 2))
        nL=size(utrL, 1)
        nR=size(utrR, 1)
        if strand[1]=='+'
            addrows!(utr5, tb(pos=utrL, chrno=fill(chrno[1], nL), strand=fill(strand[1], nL)))
            addrows!(utr3, tb(pos=utrR, chrno=fill(chrno[1], nR), strand=fill(strand[1], nR)))
        else
            @assert strand[1]=='-'
            addrows!(utr3, tb(pos=utrL, chrno=fill(chrno[1], nL), strand=fill(strand[1], nL)))
            addrows!(utr5, tb(pos=utrR, chrno=fill(chrno[1], nR), strand=fill(strand[1], nR)))
        end
    end
    utr5=value(utr5)
    utr3=value(utr3)

    # lnc=rec(exn, exn["gene_typ"].=="lincRNA")
    lnc=rec(exn, .!ismbr(exn["geneid"], cdg["geneid"]))
    # gexn=rec(exn, ismbr(exn["gene_typ"], c"protein_coding, lincRNA"))
    gexn=exn
    ge=tb(c"pos, transid"=>grpfun(x->[minimum(x[:, 1]) maximum(x[:, 2])], gexn["transid"], gexn["exn_pos"]))
    d"chrno, strand"ge=dtshift(ge["transid"], exn["transid"], d"chrno, strand"exn)
    ge["promoter_pos"]=rowfun(ge["pos"], ge["strand"]) do p, s
        if s[1]=='+'
            [p[1]-2000 p[1]-1]
        else
            @assert s[1]=='-'
            [p[2]+1 p[2]+2000]
        end
    end

    pk=jload("merged_Insm1_peak_Adult.jld2")
    # pk=jload("merged_Insm1_peak_E18.jld2")
    
    pk["group"]=fill("", rnum(pk))

    _, pi=genomemap(d"chrno, cds_pos"cdg, (pk["chrno"], pk["abs_summit"]|>x->[x x]))
    pk["group"][pi].="CDS"

    l=findall(isempty.(pk["group"]))
    _, pi=genomemap(d"chrno, pos"utr5, (pk["chrno"][l], pk["abs_summit"][l]|>x->[x x]))
    pk["group"][l[pi]].="5'-UTR"

    l=findall(isempty.(pk["group"]))
    _, pi=genomemap(d"chrno, pos"utr3, (pk["chrno"][l], pk["abs_summit"][l]|>x->[x x]))
    pk["group"][l[pi]].="3'-UTR"

    l=findall(isempty.(pk["group"]))
    _, pi=genomemap(d"chrno, promoter_pos"ge, (pk["chrno"][l], pk["abs_summit"][l]|>x->[x x]))
    pk["group"][l[pi]].="promoter"

    l=findall(isempty.(pk["group"]))
    _, pi=genomemap(d"chrno, exn_pos"lnc, (pk["chrno"][l], pk["abs_summit"][l]|>x->[x x]))
    pk["group"][l[pi]].="ncRNA_exon"

    l=findall(isempty.(pk["group"]))
    _, pi=genomemap(d"chrno, pos"ge, (pk["chrno"][l], pk["abs_summit"][l]|>x->[x x]))
    pk["group"][l[pi]].="intron"

    l=findall(isempty.(pk["group"]))
    pk["group"][l].="intergenetic"

    jsave("Insm1_peak_Adult_with_gene_structures", pk)
end

## Check
pk=jload("Insm1_peak_Adult_with_gene_structures")
freqls(pk["group"])
freqls(pk["group"][pk["smpNO"].=="12"])
#}}

#{{ Aire binding site analysis
# narrowPeak files were downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165550
fns=c"""
Aire_ChIP_WT/GSM5036536_Aire_C313Y_WT_ChIP_rep1_vs_Aire_C313Y_WT_WCE_rep1_peaks.narrowPeak.gz
Aire_ChIP_WT/GSM5036537_Aire_C313Y_WT_ChIP_rep2_vs_Aire_C313Y_WT_WCE_rep2_peaks.narrowPeak.gz
Aire_ChIP_WT/GSM5036544_C442G_WT2_ChIP_vs_C442G_WT2_WCE_peaks.narrowPeak.gz
Aire_ChIP_WT/GSM5036545_C442G_WT3_ChIP_vs_C442G_WT3_WCE_peaks.narrowPeak.gz
"""c
lbs=cut"_1".(cut"/2".(fns))
T=mapr(fns, lbs, mrows=true) do fn, lb
    T=readtb(sfc"zcat $1"(fn), head=c"chr, N::pos, <, X::, X::, X::, F::signalValue, F::pValue, F::qValue, N::peak")
    T["pos"][:, 1].+=1
    tbexp(T, GSMid=lb)
end
T["summit_pos"]=(T["pos"][:, 1].+T["peak"]) |>x->[x x]
jsave("Aire_ChIP_WT_4rep_narrowPeak", T)
#}}
#{{ Check Insm1 overlap with Aire
pk1=jload("merged_Insm1_peak_Adult_2hit.jld2")
pk2=jload("../E18_Insm1_CTCF_ATAC_analysis/merged_Insm1_peak_E18_2hit.jld2")
A=jload("../E18_Insm1_CTCF_ATAC_analysis/Aire_peak_overlap/Aire_ChIP_WT_4rep_narrowPeak")
A=rec(A, A["qValue"].>=5)

iI1, iA1=genomemap(d"chr, pos"pk1, d"chr, pos"A, touch=true)
pk1["isAire"]=falsesbut(rnum(pk1), iI1)
showprop(pk1["isAire"])

iI2, iA2=genomemap(d"chr, pos"pk2, d"chr, pos"A, touch=true)
pk2["isAire"]=falsesbut(rnum(pk2), iI2)
showprop(pk2["isAire"])

A["isInsm1"]=[falsesbut(rnum(A), iA1) falsesbut(rnum(A), iA2)]
showprop(A["isInsm1"])
#}}
