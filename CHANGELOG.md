

# Changelog

## initial release v1.0.0 2025-1-7

## RABBIT v1.1.4

* MagicImpute: additional likelihood-based accept/reject for founderimpute
* MagicImpute: new option "isallowmissing"
* MagicCall: fix a bug wehn no likeparameters to be inferred. 
* MagicMap: revise internal functions for determinning default options (of value nothing) such as isdupebinning, isrfbining, and minsilhouette
* MagicBase: update parsebreed for an input breeder's pedfile with some sample indviduals acting also as founders in pedcode
* MagicReconstruct: new option "thincm"

## RABBIT v1.2.3

* MagicFilter and others: rename options by removing "snp_". For example, change snp_minmaf to minmaf
* MagicImpute: additional likelihood-based accept/reject for foundercorrect
* MagicImpute: delay accet/reject for founderimpute, which is used only for stoping the substep of founder imputation. 
* MagicImpute: new option startbyhalf = 5; when iteration >= startbyhalf, obtain proposal for imputation for each founder block conditional on the other half chromosomes (i.e. byhalf=true)
* MagicImpute: new option threshproposal = 0.7; when obtainning proposal for imputation, founder alleles with conditional probability < threshproposal will be set to missing. 
* MagicImpute: if isallowmissing = false (default true), threshproposal is set to -1 in the last iteration and the proposal is always accepted such that all the remaining missing founder genotypes will be imputed. 
* MagicImpute: random founder partition only in the first iteration rather than in every iteration

## RABBIT v1.2.4

* magicfilter CLI:  add options minmaf and maxmiss 
* magicmap CLI: set type of knncluster and knnorder to abstractstring (parsed to function as a number of #markers)
* magiccall: revise wdict in priorfhaplo_singlesite
