

# Changelog

## initial release v1.0.0 2025-1-7

## RABBIT v1.1.4

* MagicImpute: additional likelihood-based accept/reject for founderimpute
* MagicImpute: new option "isallowmissing"
* MagicCall: fix a bug wehn no likeparameters to be inferred. 
* MagicMap: revise internal functions for determinning default options (of value nothing) such as isdupebinning, isrfbining, and minsilhouette
* MagicBase: update parsebreed for an input breeder's pedfile with some sample indviduals acting also as founders in pedcode
* MagicReconstruct: new option "thincm"

## RABBIT v1.2.1

* MagicFilter and others: rename options by removing "snp_". For example, change snp_minmaf to minmaf
* MagicImpute: additional likelihood-based accept/reject for foundercorrect
* MagicImpute: delay accet/reject for founderimpute
* MagicImpute: reset default isallowmissing = true
