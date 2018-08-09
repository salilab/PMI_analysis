import glob

files = glob.glob('data/XLs_satisfaction_all*csv')

chains = {'Vif':'A', 'CBFB':'B', 'EloB':'C', 'EloC':'D', 'A3G':'E', 'CUL5':'F', 'Rbx2':'G'}

out = open('XLs_pseudobonds.dat','w')
for f in files:
    for line in open(f,'r'):
        if 'Distance' in line:
            vals_0 = line.split(',')
            vals_1 = vals_0[0].split('|')
            if float(vals_0[1]) < 0.1:
                color = 'orange'
            elif float(vals_0[1])>0.9:
                color = 'forest green'
            elif float(vals_0[1]) >= 0.1 and float(vals_0[1])<=0.9:
                color = 'gold'

            # Get chains
            ch0 = chains[vals_1[3]]
            ch1 = chains[vals_1[5]]
        
            resi0 = vals_1[4]
            resi1 = vals_1[6]
            
            out.write("#0:"+resi0+"."+ch0+"@CA "+ "#0:"+resi1+"."+ch1+"@CA "+color+" \n")

out.close()
