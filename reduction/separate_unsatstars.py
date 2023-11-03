from saturated_star_finding import iteratively_remove_saturated_stars, remove_saturated_stars

basepath = '/orange/adamginsburg/jwst/cloudc/'
proposal_id = '2221'
field = '002'
filtername = 'F405N'
modules = ['nrca', 'nrcb', 'merged']

for module in modules: 
    try:
        remove_saturated_stars(f'{basepath}/{filtername}/pipeline/jw0{proposal_id}-o{field}_t001_nircam_clear-{filtername.lower()}-{module}_i2d.fits')
        remove_saturated_stars(f'{basepath}/{filtername}/pipeline/jw0{proposal_id}-o{field}_t001_nircam_clear-{filtername.lower()}-{module}_realigned-to-vvv.fits')
    except (TimeoutError, requests.exceptions.ReadTimeout) as ex:
        print("Failed to run remove_saturated_stars with failure {ex}")