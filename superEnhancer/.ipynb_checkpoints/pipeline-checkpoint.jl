include("../deps/dataProcessKit/dataProcessKit.jl")
using .dataProcessKit

#{{ Super enhancer positions provided by author
T=readtb("super_enhancer_mTECs_CBDM_lab.bed", head=c"chr, N::pos, <")
T["pos"][:, 1]=T["pos"][:, 1].+1
T["chrno"]=chr2no(T["chr"])
T["supEnhNO"]=collect(1:rnum(T))
jsave("super_enhancer_mTECs_providedByAuthor", T)
#}}
