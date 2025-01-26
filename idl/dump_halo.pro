; This routine dumps all particles/grids of subhalos and unbound components in target FoF halos.
; nout: snapshot number
; list: one dim array indicating the IDs of FoF halos of interest ex) [100,300,10000]


pro dump_halo,nout=nout,list=list



	h=0.684d

;	info_path='~/HR5/data/02_without_dc/info'
	out_path='catalogue/'
	snout1=string(nout,format='(i05)')
;	rd_info,info,file=info_path+'/info_'+snout1+'.txt'
;	h=info.h0/100.
	info=1.
	lbox=717.22904d0
	
	snout0=string(nout,format='(i03)')

;	path='~/HR5/psb/low_res/sn'+snout0+'/'
	path='~/HR5/psb/catalogue/sn'+snout0+'/halo_data'
	ft=file_test(path)
	if(ft[0] eq 0) then spawn,'mkdir -p '+path

	print,'Extracting data from snapshot '+snout0
	file0='FoF_Data/FoF.'+snout1+'/GALFIND.DATA.'+snout1
;	file0='FoF_Data/low_res/FoF.'+snout1+'/GALFIND.DATA.'+snout1
	file1='FoF_Data/FoF.'+snout1+'/GALCATALOG.LIST.'+snout1
	file2='FoF_Data/FoF.'+snout1+'/GALFIND.CENTER.'+snout1
	file3='FoF_Data/FoF.'+snout1+'/background_ptl.'+snout1
;	file3='FoF_Data/low_res/FoF.'+snout1+'/background_ptl.'+snout1
;	file4='catalogue/sn'+snout0+'/halo_catalogue.sav'
;	file4='FoF_Data/low_res/FoF.'+snout1+'/halo_catalogue_'+snout1+'.sav'
;	file5='catalogue/sn'+snout0+'/link_sub_host_galaxy.sav'
;	restore,file4
;	restore,file5

	thalo={nsub:0L, ndm:0L, nstar:0L, nsink: 0L, ngas: 0L, npall: 0L, mtot:0.d, mdm:0.d, mgas:0.d, msink:0.d, mstar:0.d, pos:dblarr(3), vel:dblarr(3)}
	tsub={ndm:0L, ngas:0L, nsink: 0L, nstar: 0L, npall: 0L, dum:0L, mtot:0.d, mdm:0.d, mgas:0.d, msink:0.d, mstar:0.d, pos:dblarr(3), vel:dblarr(3)}
	tpart={pos:dblarr(3),vel:dblarr(3),mass:0.d,dum0:0.d,tp:0.d,zp:0.d,mass0:0.d,tpp:0.d, indtab:0.d,id:0LL,potential:0.d,level:0L,dum1:0.}
	tsink={pos:dblarr(3),vel:dblarr(3),mass:0.d,tbirth:0.d,angm:dblarr(3),ang:dblarr(3),dmsmbh:dblarr(3),esave:0.d,smag:0.d,eps:0.d,id:0L,dum0:0L}
	tgas={pos:dblarr(3),dx:0.d,vel:fltarr(3),dum0:0.,density:0.d,temp:0.,metal:0.,fe:0.,h:0.,o:0.,level:0L,mass:0.,dum1:0.,id:0LL,potential:0.d,f:dblarr(3)}

	bhalo=n_tags(thalo,/length)
	bsub=n_tags(tsub,/length)
	bpart=n_tags(tpart,/length)
	bsink=n_tags(tsink,/length)
	bgas=n_tags(tgas,/length)
	

;===================================
; Reading halo catalogue
;===================================

	nhalo=0L
	nsub=0L
	openr,10,file1
    while ~eof(10) do begin

        nhalo=nhalo+1
        readu,10,thalo
        for i=0L,thalo.nsub-1 do begin
            nsub=nsub+1
            readu,10,tsub
        endfor
    endwhile
	close,10

	halo=replicate(thalo,nhalo)
	subhalo=replicate(tsub,nsub)

    openr,10,file
    k=0L
    for i=0L,nhalo-1 do begin
        readu,10,thalo
        halo[i]=thalo
        for j=0L,thalo.nsub-1 do begin
            readu,10,tsub
            subhalo[k]=tsub
            k=k+1
        endfor
    endfor
    close,10


	tnstar=total(halo.nstar,/integer)

	openr,10,file0
	openr,20,file3


;===============================================
;	Making tables of cumultive number of particles/grids/subhalos 
;	to compute byte locating target halos
;===============================================

	idx=where(list ge 0)
	ctnstar=total(subhalo[0:nsub-2].nstar,/integer,/cumulative)
	ctndm=total(subhalo[0:nsub-2].ndm,/integer,/cumulative)
	ctngas=total(subhalo[0:nsub-2].ngas,/integer,/cumulative)
	ctnsink=total(subhalo[0:nsub-2].nsink,/integer,/cumulative)
	ctnsub=total(halo[0:nhalo-2].nsub,/integer,/cumulative)

	ctnstar=[0,ctnstar]
	ctndm=[0,ctndm]
	ctngas=[0,ctngas]
	ctnsink=[0,ctnsink]
	ctnsub=[0,ctnsub]

;	bndm=halo.ndm
;	bnstar=halo.nstar
;	bngas=halo.ngas
;	bnsink=halo.nsink

	bndm=unbound.ndm
	bnstar=unbound.nstar
	bngas=unbound.ngas
	bnsink=unbound.nsink


	bctndm=total(bndm[0:nhalo-2],/integer,/cumulative)
	bctnstar=total(bnstar[0:nhalo-2],/integer,/cumulative)
	bctngas=total(bngas[0:nhalo-2],/integer,/cumulative)
	bctnsink=total(bnsink[0:nhalo-2],/integer,/cumulative)
	bctndm=[0,bctndm]
	bctnstar=[0,bctnstar]
	bctngas=[0,bctngas]
	bctnsink=[0,bctnsink]


;    bskip_bound=lon64arr(nhalo)
 ;   bskip_unbound=lon64arr(nhalo)
  ;  boutput='catalogue/halo_byte_skip_'+snout0+'.sav'
;	nsub_halo=halo.nsub
 ;   for i=0LL,nhalo-1 do begin
;		j=ctnsub[i]
 ;       bskip_bound[i]=i*bhalo+(ctndm[j]+ctnstar[j])*bpart+ctngas[j]*bgas+ctnsink[j]*bsink+j*bsub
  ;      bskip_unbound[i]=i*bsub+(bctndm[i]+bctnstar[i])*bpart+bctngas[i]*bgas+bctnsink[i]*bsink
  ;  endfor
;    save,file=boutput,bskip_bound,bskip_unbound,nsub_halo
;	stop

	if(idx[0] ne -1) then begin
		for i=0L,n_elements(idx)-1 do begin
			hid=list[idx[i]]
			nsub=halo[hid].nsub
			gid=ctnsub[hid]
			output=path+'/hid_'+string(hid,format='(i07)')+'.sav'
			tnstar=ctnstar[gid]
			tndm=ctndm[gid]
			tngas=ctngas[gid]
			tnsink=ctnsink[gid]

			btndm=bctndm[hid]
			btnstar=bctnstar[hid]
			btnsink=bctnsink[hid]
			btngas=bctngas[hid]

			skip=long64(hid+1)*bhalo+long64(gid)*bsub+(tnstar+tndm)*bpart+tngas*bgas+tnsink*bsink
			skip1=long64(hid+1)*bsub+(btnstar+btndm)*bpart+btngas*bgas+btnsink*bsink
			print,'Halo ID : '+string(hid,format='(i07)')+' Skip bytes : '+string(skip,format='(i15)')
			print,' Skip bytes : '+string(skip1,format='(i15)')
			adm=-1
			astar=-1
			agas=-1
			asink=-1
			if(halo[hid].ndm gt 0) then adm=replicate(tpart,halo[hid].ndm)
			if(halo[hid].nstar gt 0) then astar=replicate(tpart,halo[hid].nstar)
			if(halo[hid].ngas gt 0) then agas=replicate(tgas,halo[hid].ngas)
			if(halo[hid].nsink gt 0) then asink=replicate(tsink,halo[hid].nsink)
			idm=lonarr(2,nsub+1)
			istar=lonarr(2,nsub+1)
			igas=lonarr(2,nsub+1)
			isink=lonarr(2,nsub+1)
			isub=lonarr(nsub+1)
			isub[0]=-1
			idm[*]=-1
			istar[*]=-1
			isink[*]=-1
			igas[*]=-1

			if(unbound[hid].ndm gt 0) then begin
				idm[0,0]=0
				idm[1,0]=unbound[hid].ndm-1
			endif
			if(unbound[hid].nstar gt 0) then begin
				istar[0,0]=0
				istar[1,0]=unbound[hid].nstar-1
			endif
			if(unbound[hid].ngas gt 0) then begin
				igas[0,0]=0
				igas[1,0]=unbound[hid].ngas-1
			endif
			if(unbound[hid].nsink gt 0) then begin
				isink[0,0]=0
				isink[1,0]=unbound[hid].nsink-1
			endif
			lidm=0L
			listar=0L
			ligas=0L
			lisink=0L

			point_lun,20,skip1
			if(unbound[hid].ndm gt 0) then begin
				dm=replicate(tpart,unbound[hid].ndm)
				readu,20,dm
				adm[0:idm[1,0]]=dm
				lidm=idm[1,0]+1
			endif
			if(unbound[hid].ngas gt 0) then begin
				gas=replicate(tgas,unbound[hid].ngas)
				readu,20,gas
				agas[0:igas[1,0]]=gas
				ligas=igas[1,0]+1
			endif
			if(unbound[hid].nsink gt 0) then begin
				sink=replicate(tsink,unbound[hid].nsink)
				readu,20,sink
				asink[0:isink[1,0]]=sink
				lisink=isink[1,0]+1
			endif
			if(unbound[hid].nstar gt 0) then begin
				star=replicate(tpart,unbound[hid].nstar)
				readu,20,star
				astar[0:istar[1,0]]=star
				listar=istar[1,0]+1
			endif

	
			point_lun,10,skip
			for j=0,nsub-1 do begin
				isub[j+1]=gid
				readu,10,tsub
		;		skip=skip+bsub
				if(subhalo[gid].ndm gt 0) then begin
					dm=replicate(tpart,subhalo[gid].ndm)
					readu,10,dm
					adm[lidm:lidm+subhalo[gid].ndm-1]=dm
					idm[0,j+1]=lidm
					lidm=lidm+subhalo[gid].ndm
					idm[1,j+1]=lidm-1
				endif

				if(subhalo[gid].ngas gt 0) then begin
					gas=replicate(tgas,subhalo[gid].ngas)
					readu,10,gas
					gas.temp=gas.temp/gas.density
					agas[ligas:ligas+subhalo[gid].ngas-1]=gas
					igas[0,j+1]=ligas
					ligas=ligas+subhalo[gid].ngas
					igas[1,j+1]=ligas-1
				endif	

				if(subhalo[gid].nsink gt 0) then begin
					sink=replicate(tsink,subhalo[gid].nsink)
					readu,10,sink
					asink[lisink:lisink+subhalo[gid].nsink-1]=sink
					isink[0,j+1]=lisink
					lisink=lisink+subhalo[gid].nsink
					isink[1,j+1]=lisink-1
				endif
				if(subhalo[gid].nstar gt 0) then begin
					star=replicate(tpart,subhalo[gid].nstar)
					readu,10,star
					astar[listar:listar+subhalo[gid].nstar-1]=star
					istar[0,j+1]=listar
					listar=listar+subhalo[gid].nstar
					istar[1,j+1]=listar-1
				endif
;				skip_lun,10,bsub
				gid=gid+1
			endfor

			save,file=output,adm,astar,agas,asink,info,hid,idm,igas,isink,istar,isub

		endfor
	endif
	
	
	close,10
	close,20
end
