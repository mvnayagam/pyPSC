from copy import Error
import h5py
import re
import numpy as np
import polytope as pc


def wrtdata_v0(fname, rc, volume, err, final, extremepnts,volAsym, Lsol):
    
    with h5py.File(fname, 'a') as f:
    ### Writing volume information in file 
        if 'lenofsolution' in f:
            lsolu = str('/lenofsolution/')+str('lsol') +str(rc)
            f.create_dataset(lsolu,  data=Lsol,  dtype='float64')
        else:
            gvol = f.create_group('lenofsolution')
            lsolu = str('/lenofsolution/')+str('lsol') +str(rc)
            f.create_dataset(lsolu,  data=Lsol,  dtype='float64')
    
    ### Writing volume information in file 
        if 'vol' in f:
            v = str('/vol/')+str('v') +str(rc)
            f.create_dataset(v,  data=volume,  dtype='float64')
        else:
            gvol = f.create_group('vol')
            v = str('/vol/')+str('v') +str(rc)
            f.create_dataset(v,  data=volume,  dtype='float64')
    
    ### Writing error information in file 
        if 'error' in f:
            derr = str('/error/')+str('err') +str(rc)
            f.create_dataset(derr,  data=err,  dtype='float64')
        else:
            gerror = f.create_group('error')
            derr = str('/error/')+str('err') +str(rc)
            f.create_dataset(derr,  data=err,  dtype='float64')
    
    ### Writing polytope information in file 
        if 'polytope' in f:
            pa=str('/polytope/')+str('pa') +str(rc)
            pb=str('/polytope/')+str('pb') +str(rc)
            f.create_dataset(pa, data=final.A, dtype='float64')
            f.create_dataset(pb, data=final.b, dtype='float64')
        else:
            gpolytope = f.create_group('polytope')
            pa=str('/polytope/')+str('pa') +str(rc)
            pb=str('/polytope/')+str('pb') +str(rc)
            f.create_dataset(pa, data=final.A, dtype='float64')
            f.create_dataset(pb, data=final.b, dtype='float64')
    
    ### Writing extreme points of polytope information in file 
        if 'extreme' in f:
            ext=str('/extreme/')+str('pa') +str(rc)
            f.create_dataset(ext, data=extremepnts, dtype='float64')
        else:
            gpolytope = f.create_group('extreme')
            ext=str('/extreme/')+str('pa') +str(rc)
            f.create_dataset(ext, data=extremepnts, dtype='float64')
    
    ### Writing sum of all volums in Asym part 
        if 'total_volume_in_Asym' in f:
            v_asym = str('/total_volume_in_Asym/')+str('pair') +str(rc)
            f.create_dataset(v_asym,  data=volAsym,  dtype='float64')
        else:
            gerror = f.create_group('total_volume_in_Asym')
            v_asym = str('/total_volume_in_Asym/')+str('pair') +str(rc)
            f.create_dataset(v_asym,  data=volAsym,  dtype='float64')        
    return

def wrtcoor(fname, pairs):
    
    ### Writing polygenerated coordinats to file
        
    with h5py.File(fname, 'a') as f:
        for i, Pair in enumerate(pairs):
            Pair_sort = np.sort(Pair)[::-1]
            
            if 'generatedcoordinate' in f:
                co = str('/generatedcoordinate/')+str(i)
                f.create_dataset(co, data=Pair,  dtype='float64')
            else:
                gco = f.create_group('generatedcoordinate')
                co  = str('/generatedcoordinate/')+str(i)
                f.create_dataset(co, data=Pair,  dtype='float64')
            
            # if 'coordinate' in f:
            #     co = str('/coordinate/')+str(i)
            #     f.create_dataset(co, data=Pair,  dtype='float64')
            # else:
            #     gco = f.create_group('coordinate')
            #     co  = str('/coordinate/')+str(i)
            #     f.create_dataset(co, data=Pair,  dtype='float64')
            
            if 'coordinate_sorted' in f:
                co = str('/coordinate_sorted/')+str(i)
                f.create_dataset(co, data=Pair_sort, dtype='float64')
            else:
                gco = f.create_group('coordinate_sorted')
                co  = str('/coordinate_sorted/')+str(i)
                f.create_dataset(co, data=Pair_sort, dtype='float64')
    return

def wrtvolume(fname, rc, volume, dx, dy, final):
    
    with h5py.File(fname, 'a') as f:
        v =str('v') +str(rc)
        dxx=str('dx')+str(rc)
        dyy=str('dy')+str(rc)
        
        f.create_dataset(v,  data=volume,  dtype='float64')
        f.create_dataset(dxx, data=dx, dtype='float64')
        f.create_dataset(dyy, data=dy, dtype='float64')
        
        pa=str('s')+str(rc)+str('a')
        pb=str('s')+str(rc)+str('b')
        
        f.create_dataset(pa, data=final.A, dtype='float64')
        f.create_dataset(pb, data=final.b, dtype='float64')
    return

def wrtallsolution(fname, rc, solutionall):
        
    def wrtfile(rc, file, sa):
        count=0
        for cou, ip in enumerate(sa):
            if type(ip) is pc.Polytope:
                
                allpa=str('/allsolution_polytope/Pair')+ str(rc)+ str('/a')+ str(count)
                allpb=str('/allsolution_polytope/Pair')+ str(rc)+ str('/b')+ str(count)
                
                file.create_dataset(allpa, data=ip.A, dtype='float64')
                file.create_dataset(allpb, data=ip.b, dtype='float64')
                count += 1
                
            elif type(ip) is pc.Region:
                for iqinx, iq in enumerate(ip):
                    count += 1
                    #cou = iqinx + cou
                    allpa=str('/allsolution_polytope/')+ str('/Pair/')+ str(rc)+ str('/a')+ str(rc)+ str(count)
                    allpb=str('/allsolution_polytope/')+ str('/Pair/')+ str(rc)+ str('/b')+ str(rc)+ str(count)
                    file.create_dataset(allpa, data=iq.A, dtype='float64')
                    file.create_dataset(allpb, data=iq.b, dtype='float64')
        
    with h5py.File(fname, 'a') as file:
        
        #---> Writing extreme points of polytope information in file 
        if 'allsolution_polytope' in file:
            grp=file['/allsolution_polytope/']
            sgp=str('Pair')+ str(rc)
            sg = grp.create_group(sgp)
            wrtfile(rc, file, solutionall)
            
        else:
            gpolytope = file.create_group('allsolution_polytope')
            spg=str('Pair')+ str(rc)
            sg = gpolytope.create_group(spg)
            
            wrtfile(rc, file, solutionall)
        
    return

def wrttime_mc(rc, fname, timeinfo):
    with h5py.File(fname, 'a') as f:
        if 'time_total' in f:
            dtime = str('/time_total/')+str('tforPair') +str(rc)
            f.create_dataset(dtime,  data=timeinfo,  dtype='float64')
        else:
            gtime = f.create_group('time_total')
            dtime = str('/time_total/')+str('tforPair') +str(rc)
            f.create_dataset(dtime,  data=timeinfo,  dtype='float64')
        
    return

def readoldsolution_v0(pairID, fname):
    
    with h5py.File(fname, 'r') as file:
        
        p=file.get('allsolution/Pair'+str(pairID))
        polyold=[]
        
        for ic, ds in enumerate(p.keys()):
            if ic < int(len(p.keys())/2):
                db='b'+re.split('(\d+)',ds)[1]
                da=np.array(p[ds])
                polyold.append(pc.Polytope(np.array(p[ds]), np.array(p[db])))
    
    return polyold

def readh5file(fn):
    
    with h5py.File(fn, 'r') as f:
        mean = []
        ls = list(f.items())
        
        vols = f.get('vol')
        vall = [np.array(vols.get(i)) for i in np.array(vols)]
        mean.append(np.mean(np.array(vall)))
        
        er=f.get('error')
        error=[np.array(er.get(i)) for i in np.array(er)]
                
        koorkey  = f.get('generatedcoordinate')
        koor = [np.array(koorkey.get(i)) for i in np.array(koorkey)]
        
        unsortkoorkey=f.get('coordinate')
        unsortkoor=[np.array(unsortkoorkey.get(i)) for i in np.array(unsortkoorkey)]
        
        ext=f.get('extreme')
        extre=[np.mean(ext.get(i),0) for i in np.array(ext)]
        
        radi=f.get('total_volume_in_Asym')
        radius=[np.array(radi.get(i)[1]) for i in np.array(radi)]
        
        
        noofsolu=f.get('allsolution')
        noofsolution = [ (np.shape( np.array(f.get('allsolution/Pair'+str(ii))) )[0])/2 for ii in range(1, len(noofsolu)+1) ]
        
    return np.array(vall), mean[0], error, np.array(koor),  np.array(unsortkoor), extre, radius, noofsolution

def readh5file_v2(fn):
    
    with h5py.File(fn, 'r') as f:
        mean = []
        ls = list(f.items())
        
        vol = f.get('solution_volume')
        volume = [np.array(vol.get(i)) for i in np.array(vol)]
        mean.append(np.mean(np.array(volume)))
        
        koorkey  = f.get('generatedcoordinate')
        koordinate = [np.array(koorkey.get(i)) for i in np.array(koorkey)]
        
        err  = f.get('solution_error')
        error = [np.array(err.get(i)) for i in np.array(err)]
                
        ext=f.get('solution_extremepoint')
        centroid=[np.mean(ext.get(i),0) for i in np.array(ext)]
        
        radkey=f.get('grandradius')
        radius=[np.array(radkey.get(i)) for i in np.array(radkey)]
        
        noofsolu=f.get('allsolution')
        noofsolution = [ (np.shape( np.array(f.get('allsolution/Pair'+str(ii))) )[0])/1 for ii in range(0, len(noofsolu)) ]
        
        noofsolu=f.get('solution_totalNr')
        noofsolution=[np.array(noofsolu.get(i)) for i in np.array(noofsolu)]
        
    return np.array(volume), mean[0], error, mean, np.array(koordinate), centroid, radius, noofsolution

def readoldsolution(pairID, fname):
    
    with h5py.File(fname, 'r') as file:
        
        p=file.get('allsolution_polytope/Pair'+str(pairID))
        polyold=[]
        
        for ic, ds in enumerate(p.keys()):
            if ic < int(len(p.keys())/2):
                db='b'+re.split('(\d+)',ds)[1]
                da=np.array(p[ds])
                polyold.append(pc.Polytope(np.array(p[ds]), np.array(p[db])))
    
    return polyold

def wrtdata(pair_inx, fname, solution, solution_polytope, solution_volume, solution_error, solution_extremepnts, vol_in_Asym, grandradius , total_solNr, allsolutions):
    
    with h5py.File(fname, 'a') as f:
        
        # ---> Writing found solution for given structure in file
        if 'solution' in f:
            solkey = str('/solution/')+str('solutionforPair') +str(pair_inx)
            f.create_dataset(solkey, data=solution, dtype='float64')
        else:
            gvol   = f.create_group('solution')
            solkey = str('/solution/')+str('solutionforPair') +str(pair_inx)
            f.create_dataset(solkey, data=solution, dtype='float64')

        # ---> Writing all found solutions in file
        if 'allsolution' in f:
            allsolution = str('/allsolution/')+str('Pair') +str(pair_inx)
            f.create_dataset(allsolution,  data=allsolutions,  dtype='float64')
        else:
            gerror = f.create_group('allsolution')
            allsolution = str('/allsolution/')+str('Pair') +str(pair_inx)
            f.create_dataset(allsolution, data=allsolutions, dtype='float64')   
        
        # ---> Writing polytope contains given structure in file 
        if 'solution_polytope' in f:
            pa=str('/solution_polytope/')+str('pa') +str(pair_inx)
            pb=str('/solution_polytope/')+str('pb') +str(pair_inx)
            f.create_dataset(pa, data=solution_polytope.A, dtype='float64')
            f.create_dataset(pb, data=solution_polytope.b, dtype='float64')
        else:
            gpolytope = f.create_group('solution_polytope')
            pa=str('/solution_polytope/')+str('pa') +str(pair_inx)
            pb=str('/solution_polytope/')+str('pb') +str(pair_inx)
            f.create_dataset(pa, data=solution_polytope.A, dtype='float64')
            f.create_dataset(pb, data=solution_polytope.b, dtype='float64')

        # ---> Writing volume of polytope contains given structure in file
        if solution_volume:
            if 'solution_volume' in f:
                v = str('/solution_volume/')+str('v') +str(pair_inx)
                f.create_dataset(v,  data=solution_volume,  dtype='float64')
            else:
                gvol = f.create_group('solution_volume')
                v = str('/solution_volume/')+str('v') +str(pair_inx)
                f.create_dataset(v,  data=solution_volume,  dtype='float64')

        # ---> Writing error information in file
        if np.all(solution_error):
            if 'solution_error' in f:
                derr = str('/solution_error/')+str('err') +str(pair_inx)
                f.create_dataset(derr,  data=solution_error,  dtype='float64')
            else:
                gerror = f.create_group('solution_error')
                derr = str('/solution_error/')+str('err') +str(pair_inx)
                f.create_dataset(derr,  data=solution_error,  dtype='float64')
        
        # ---> Writing extreme points of polytope contains given structure in file 
        if 'solution_extremepoint' in f:
            ext=str('/solution_extremepoint/')+str('pa') +str(pair_inx)
            f.create_dataset(ext, data=solution_extremepnts, dtype='float64')
        else:
            gpolytope = f.create_group('solution_extremepoint')
            ext=str('/solution_extremepoint/')+str('pa') +str(pair_inx)
            f.create_dataset(ext, data=solution_extremepnts, dtype='float64')
                
        # ---> Writing sum of all volume in asymmetric PS 
        if vol_in_Asym:
            if 'total_volume_in_Asym' in f:
                v_asym = str('/total_volume_in_Asym/')+str('pair') +str(pair_inx)
                f.create_dataset(v_asym,  data=vol_in_Asym,  dtype='float64')
            else:
                gerror = f.create_group('total_volume_in_Asym')
                v_asym = str('/total_volume_in_Asym/')+str('pair') +str(pair_inx)
                f.create_dataset(v_asym,  data=vol_in_Asym,  dtype='float64')  
        
        # ---> Writing grand radius information in file
        if grandradius:
            if 'grandradius' in f:
                radius = str('/grandradius/')+str('r') +str(pair_inx)
                f.create_dataset(radius,  data=grandradius,  dtype='float64')
            else:
                gradius = f.create_group('grandradius')
                radius = str('/grandradius/')+str('r') +str(pair_inx)
                f.create_dataset(radius,  data=grandradius,  dtype='float64')
                
        # ---> Writing total no of solution information in file
        if  total_solNr:
            if 'solution_totalNr' in f:
                lsolu = str('/solution_totalNr/')+str('lsol') +str(pair_inx)
                f.create_dataset(lsolu,  data=total_solNr,  dtype='float64')
            else:
                gsol  = f.create_group('solution_totalNr')
                lsolu = str('/solution_totalNr/')+str('lsol') +str(pair_inx)
                f.create_dataset(lsolu,  data=total_solNr,  dtype='float64')
    
    return
