for filter in F410M F405N F466N; do
    for module in nrca nrcb; do
        for dao in "--daophot --skip-crowdsource" " "; do
            sbatch --array=0-23 --job-name=webb-cat-${filter}-${module}-eachexp-cloudc --output=webb-cat-${filter}-${module}-eachexp-cloudc_%j-%A_%a.log  --account=adamginsburg --qos=adamginsburg-b --ntasks=2 --nodes=1 --mem=16gb --time=96:00:00 --wrap "/blue/adamginsburg/adamginsburg/miniconda3/envs/python310/bin/python /orange/adamginsburg/jwst/cloudc/analysis/crowdsource_catalogs_long.py --filternames=${filter} --modules=${module} --each-exposure ${dao} --target=cloudc --each-suffix=destreak_o002_crf"
        done
    done
done

for filter in F212N F182M F187N; do
    for modnum in 1 2 3 4; do
        for module in nrca${modnum} nrcb${modnum}; do
            for dao in "--daophot --skip-crowdsource" " "; do
                sbatch --array=0-23 --job-name=webb-cat-${filter}-${module}-eachexp-cloudc --output=webb-cat-${filter}-${module}-eachexp-cloudc_%j-%A_%a.log  --account=adamginsburg --qos=adamginsburg-b --ntasks=2 --nodes=1 --mem=16gb --time=96:00:00 --wrap "/blue/adamginsburg/adamginsburg/miniconda3/envs/python310/bin/python /orange/adamginsburg/jwst/cloudc/analysis/crowdsource_catalogs_long.py --filternames=${filter} --modules=${module} --each-exposure ${dao} --target=cloudc --each-suffix=destreak_o002_crf"
            done
        done
    done
done