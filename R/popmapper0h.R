# file popmapper/R/popmapper0.R
# copyright (C) 2022 Hiroshi C. Ito
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License Version 2 as 
# published by the Free Software Foundation.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# A copy of the GNU General Public License is available at
# http://www.r-project.org/Licenses/

#' Experimental R package.
#' @aliases  popmapper popmapper-package
#' #' @importFrom foreach %do%
#' @keywords internal
"_PACKAGE"

flag_envs=FALSE;

##' @author Hiroshi C. Ito
#' @export
param0=list(
    flag_cols_cerwave=0,
    flag_fix_kdis=1,
    flag_welch_t_test=0,
    kdis0=2,
    flag_random=1,
    flag_perm=1, ## 1 or -1 (flip former and latter)
    nperm=5000,
    edge2ampsm=10.0,
    ngroup_show=-1, ## -1: optimize
    edge_p_value=0.001,
    flag_fix_lim=0,
    kdis0=2,
    ndiv=100,
    group_max=5,
    flag_remove_edge_null=0,
    edge_null=1000,
    legend_side="topleft",
    flag_ind_wise=1,
    flag_read_sample_location=1,
    flag_show_legend_location=1,
    flag_runid=1,
    flag_binary=2, ## 0 euclid, 1 my binary, 2 hybrid, 3 binary by R function
    nsamp=200,
    win_width=c(6,6,6,6),
    win_height=c(5,5,5,5),
    flag_pforeach=0
    );



##win_width=6;
##win_height=5;
##win_width=5;
##win_height=4;

##' @title xinit". 
#' @export
xinit <- function(win_width=param0$win_width,win_height=param0$win_height){
##    win_width=param$win_width;
 ##   win_height=param$win_height;
    
    graphics.off();
    
    X11(width=win_width[1],height=win_height[1]);
    X11(width=win_width[2],height=win_height[2]);
    X11(width=win_width[3],height=win_height[3]);
    X11(width=win_width[4],height=win_height[4]);
}

fn <- function(x,y,mx,my,sx,sy){
     v=(1.0/(2*pi*sx*sy))*exp(-(x-mx)^2/(2*sx^2))*exp(-(y-my)^2/(2*sy^2));
     return(v);
 }


#' @export
find_peak <- function(pcoa,ampsm,ndiv=100){
    px=pcoa$points[,1]
    py=pcoa$points[,2]
##    ampsm=0.08;
    sx=ampsm*(max(px)-min(px));
    sy=ampsm*(max(py)-min(py));
    ##ndiv=80;
    xx0=seq(min(px)-sx,max(px)+sx,,ndiv);
    yy0=seq(min(py)-sy,max(py)+sy,,ndiv);
    xx=matrix(xx0, nrow=ndiv, ncol=ndiv, byrow=F); 
    yy=matrix(yy0, nrow=ndiv, ncol=ndiv, byrow=T); 
    
    zz=yy*0;
    
    zz=zz*0.0;
    for(k in 1:length(px)){
        zz=zz+fn(xx,yy,px[k],py[k],sx,sy);
    }

    difmask=c(1,seq(ndiv),ndiv);
    ndiv=length(zz[,1]);
    wx0=(zz>zz[difmask[1:ndiv],]);
    wx1=(zz>zz[difmask[3:(ndiv+2)],]);
    wy0=(zz>zz[,difmask[1:ndiv]]);
    wy1=(zz>zz[,difmask[3:(ndiv+2)]]);
    
    wxy0=(zz>zz[difmask[1:ndiv],difmask[3:(ndiv+2)]]);
    wxy1=(zz>zz[difmask[3:(ndiv+2)],difmask[1:ndiv]]);
    wxy2=(zz>zz[difmask[1:ndiv],difmask[1:ndiv]]);
    wxy3=(zz>zz[difmask[3:(ndiv+2)],difmask[3:(ndiv+2)]]);
    
    maskp=(wx0+wx1+wy0+wy1+wxy0+wxy1+wxy2+wxy3==8);
    lisp=which(maskp>0);
    peak=t(which(maskp>0,arr.ind=TRUE));
    peakx=xx[lisp];
    peaky=yy[lisp];
    return(list(peak=peak,peakx=peakx,peaky=peaky,lisp=lisp,ngroup=length(lisp),xx0=xx0,yy0=yy0,xx=xx,yy=yy,zz=zz,sx=sx,sy=sy));
}


#' @export
calc_dist_unit <- function(i,nn,tab11,mask,amp_euc,amp_bin,flag_binary=2){
        ##     for(i in 1:nn){
        if(i%%50==0)cat(100*i/nn,"% ");
        idm00=rep(i,nn);
        idm11=seq(nn);
        mask=tab11[i,]>0;
        nm=sum(mask);
        if(flag_binary==1){
            vv=(tab11[idm00,mask]>0)-(tab11[idm11,mask]>0);
            vv=rowSums(abs(vv)>0)/nm;
        }
        if(flag_binary==0){
            vv=(tab11[idm00,mask])-(tab11[idm11,mask]);
            vv=rowSums(abs(vv)>0)/nm;
        }
        if(flag_binary==2){
            vv00=tab11[idm00,mask];
            vv11=tab11[idm11,mask];
            vv_bi=(vv00>0)-(vv11>0);
            vv_eu=(vv00)-(vv11);
            vv_bi=(abs(vv_bi)>0);
            vv_eu=(abs(vv_eu)>0);
            vv_eu=vv_eu-vv_bi;
            n_eu_eff=rowSums((vv00>0)+(vv11>0)==2);
            vv=amp_euc*rowSums(vv_eu)/nm+amp_bin*rowSums(vv_bi)/nm;
        }
            return(c(vv,mean(n_eu_eff)));
}


calc_dist <- function(tab1,amp_euc=1,amp_bin=0,flag_binary=2,flag_pforeach=1){
    nn=nrow(tab1);
    nc=ncol(tab1);
    
    tab11=as.matrix(tab1);
    
    p_NULL=rowSums(tab11==0)/nc;
    p_NULL_mean=mean(p_NULL);
    cat("fraction of NULLs   mean:",p_NULL_mean,"   max:",max(p_NULL),"   min:",min(p_NULL), "\n");
    
    x=seq(-1,1,,nn); 
    dist0=matrix(rep(0,nn), nrow=nn, ncol=nn, byrow=T);
    dist0eu=dist0;
    dist0bi=dist0;
    
    idm0=dist0;
    idm1=dist0;
    
    for(i in 1:nn){
        for(j in 1:nn){
            idm0[i,j]=i;
            idm1[i,j]=j;
        }    
    }
    
    idm0=as.numeric(idm0);
    idm1=as.numeric(idm1);
    
    
    nn=nrow(tab11);

    if(flag_pforeach==1){
        print("using pforeach..");
        dist0=pforeach::pforeach(i=1:nn, .combine='rbind',.errorhandling="stop")({    
            calc_dist_unit(i,nn,tab11,mask,amp_euc,amp_bin,flag_binary=flag_binary);
        });
    }else{
        print("using foreach..");
        dist0=foreach::foreach(i=1:nn,.combine='rbind',.packages="foreach")%do%{    
            calc_dist_unit(i,nn,tab11,mask,amp_euc,amp_bin,flag_binary=flag_binary);
        };
    }

    
    n_eu=dist0[,ncol(dist0)];
    dist0=dist0[,-ncol(dist0)];
    dist0=0.5*(dist0+t(dist0));
    return(list(neu=n_eu,dist=dist0));
}

        
sum_ind_wise <- function(dist0,idsu){
    nnu=length(idsu);
    
    dist0u=matrix(rep(0,nnu), nrow=nnu, ncol=nnu, byrow=T);     
    for(i in 1:nnu){
        for(j in 1:nnu){
            i0=2*i-1;
            i1=2*i;
            j0=2*j-1;
            j1=2*j;            
            dist0u[i,j]=min(dist0[c(i0,i1),c(j0,j1)]);
            ##dist0u[i,j]=mean(dist0[c(i0,i1),c(j0,j1)]);
        }
    }
    return(dist0u);
}







  
read_data <- function(file_tsv=file_tsv,file_sample_location=file_sample_location,omit_popid=omit_popid,amp_euc=1,amp_bin=0,param=param0){
    if(flag_envs).ee.append("read_data",environment());
    
    print("reading file...");    
    df=read.csv(file_tsv,header=F,sep="\t",comment.char="#");
        
    tab0=df[-1,]
    tab0=as.matrix(tab0[,-c(1,2)]);
    tab0=apply(tab0,2,as.numeric);
    
    if(param$flag_random==1)tab0=tab0[,order(rnorm(ncol(tab0)))];
        
    inds=df[,1];
    inds=as.character(inds[-1]);
    sample_name=inds;
    
    if(param$flag_read_sample_location==1){
        sampleinfo=read.table(file_sample_location,header=TRUE,sep="\t");
        popids=sampleinfo$PopID;
        runids=sampleinfo$RunID;
        runlabels=as.character(sampleinfo$RunLabel);
        popnames=as.character(sampleinfo$PopName);
        
        ids1=as.integer(gsub("S","",sampleinfo$SampleName))+as.integer(sampleinfo$RunID)*1000;

        runid=unique(runids);
        rids=NULL;
        for(i in 1:length(runid)){
            rids=c(rids,(which(runids==runid[i]))[1]);
        }
        
        runlabel=runlabels[rids];
        runnames=runlabel;
        
        runnames1=paste0(runnames,"S");
        ##runid_names=paste0("S",seq(length(runnames)));
        runid_names=paste0("S",runid);
        
        for(k in 1:length(runnames1))inds=(gsub(runnames1[k],runid_names[k],inds));
        ids=as.integer(gsub("S","",inds));
        
        ##if(flag_single_run==1)ids=ids+1000;
        if(max(ids)<1000)ids=ids+1000;
        
        
        ##SampleName2=paste0(runnames[sampleinfo$RunID],sampleinfo$SampleName);
        SampleName2=paste0(sampleinfo$RunLabel,sampleinfo$SampleName);
        sampleinfo=data.frame(SampleName2=SampleName2,sampleinfo);
        
    }else{
        ids1=unique(ids);
        popids=ids1;
        runids=rep(1,length(ids1));
    }
    
    


    popids1=rep(0,nrow(tab0));
    runids1=popids1;
    mask0=rep(FALSE,nrow(tab0));
    
    for(i in 1:length(mask0)){
        lis=which(ids1==ids[i]);
        if(length(lis)>0){
            mask0[i]=TRUE;
            popids1[i]=popids[lis];
            runids1[i]=runids[lis];
            }
        }
        
    tab1=tab0[mask0,];
    ids=ids[mask0];
    popids1=popids1[mask0];
    runids1=runids1[mask0];

    
    if(length(omit_popid)>0){
        omit_mask_popid=rep(FALSE,length(popids1));
        for(j in 1:length(omit_popid))omit_mask_popid=omit_mask_popid+(popids1==omit_popid[j]);
        
        omit_mask=(omit_mask_popid>0);
        
        ids_omit=ids[omit_mask];
        popids1_omit=popids1[omit_mask];
        runids1_omit=runids1[omit_mask];
        sample_omit=sample_name[omit_mask];
        
        cat("omit_popid: ",omit_popid,"\n");
        
        cat("omitted files:\n");
        print(paste(unique(sample_omit),collapse=", "));
        
        tab1=tab0[!omit_mask,];
        ids=ids[!omit_mask];
        popids1=popids1[!omit_mask];
        runids1=runids1[!omit_mask];
        sample_name=sample_name[!omit_mask];

    }

    omit=rep(0,nrow(sampleinfo));
    if(length(omit_popid)>0){
        for(j in 1:length(omit_popid))omit=omit+(popids==omit_popid[j]);
    }
    
    sampleinfo=data.frame(sampleinfo,omit=omit);

    
    if(param$flag_remove_edge_null==1){
        mask1=(colSums(tab1==0)/nrow(tab1)<param$edge_null);    
        tab1=tab1[,mask1];
    }

    tab00=tab1;
    
    nc=ncol(tab00);
    nst=nc/2;
    if(param$flag_perm==1){
        tab1p=tab00[,(nst+1):ncol(tab1)];
        tab1=tab00[,1:nst];        
    }
    if(param$flag_perm==-1){       
        tab1p=tab00[,1:nst];
        tab1=tab00[,(nst+1):ncol(tab1)];                
    }

    return(list(tab00=tab00,tab1=tab1,tab1p=tab1p,ids=ids,popids1=popids1,runids1=runids1,popids=popids,popnames=popnames,runnames=runnames,runids=runids,runid=runid,sample.name=sample_name,sample.info=sampleinfo));
}


##' @title popmap.readdata". 
#' @export
popmap.readdata <- function(file_tsv=NULL,file_sample_location=NULL,  sample_group_name="Samples", omit_popid=NULL,amp_euc=1,amp_bin=0,param=param0){
    if(flag_envs).ee.append("popmap.readdata",environment());

  
    
    re_read_data=read_data(file_tsv=file_tsv,file_sample_location=file_sample_location, omit_popid=omit_popid,amp_euc=amp_euc,amp_bin=amp_bin,param=param0);

        tab00=re_read_data$tab00;
        tab1=re_read_data$tab1;
        tab1p=re_read_data$tab1p;
        ids=re_read_data$ids;
        popids1=re_read_data$popids1;
        popids=re_read_data$popids;
    popnames=re_read_data$popnames;
    runnames=re_read_data$runnames;
    runid=re_read_data$runid;


    
    sample.name=re_read_data$sample.name;
    runids1=re_read_data$runids1;
        ##rm(re_read_data);
        
        
        print("calculating distance...");
        
        
        res_calc_dist=calc_dist(tab1,amp_euc=amp_euc,amp_bin=amp_bin,flag_binary=param$flag_binary,flag_pforeach=param$flag_pforeach);
        res_calc_distp=calc_dist(tab1p,amp_euc=amp_euc,amp_bin=amp_bin,flag_binary=param$flag_binary,flag_pforeach=param$flag_pforeach);
        
        dist0=res_calc_dist$dist;
        n_eu=res_calc_dist$neu;    
        dist0p=res_calc_distp$dist;
        n_eup=res_calc_distp$neu;
        
        cat("\n Average loci number of both non-NULL between two samples:",mean(n_eu), "\n");
        
        dist2=dist0;
        ids2=ids;
        popids2=popids1;
        
        if(param$flag_ind_wise==1){   
            idsu=unique(ids);    
            dist0u=sum_ind_wise(dist0,idsu);
            dist2=dist0u;
            ids2=idsu;
            popids2=popids1[seq(2,length(popids1),2)];
            if(abs(param$flag_perm)==1){
                dist2p=sum_ind_wise(dist0p,idsu);        
            }
        }


    ##nrun=as.integer(max(ids2)/1000);
    nrun=max(runid);
    masks=matrix(ids2 > -1, nrow=nrun, ncol=length(ids2), byrow=T);
    for(i in 1:nrun)masks[i,]=(i*1000<=ids2)+(ids2<(i+1)*1000)==2;
    mask=colSums(masks)>0;

    
    
        popid_u=unique(popids);
        popname_u=rep("",length(popid_u));
        for(i in 1:length(popid_u))popname_u[i]=popnames[(which(popids==popid_u[i]))[1]];
        
        popnames2=popname_u[popids2];

        sample.name2=rep("a",length(ids2));
        for(i in 1:length(ids2)){
            lis=which(ids==ids2[i]);
            sample.name2[i]=sample.name[lis[1]];
        }

        
        return(list(dist2=dist2,dist2p=dist2p,ids2=ids2,popids2=popids2,popnames2=popnames2,sample.name2=sample.name2,sample.name=sample.name,popid.legend=popid_u,popname.legend=popname_u,tab00=tab00,tab1=tab1,tab1p=tab1p,ids=ids,popids1=popids1,runids1=runids1,popids=popids,popnames=popnames,neu=n_eu,neup=n_eup,masks=masks,mask=mask,nrun=nrun, runnames=runnames,runid=runid,sample.info=re_read_data$sample.info,amp_euc=amp_euc,amp_bin=amp_bin,omit_popid=omit_popid,sample_group_name=sample_group_name));   
}

get_pcoa2d <- function(dist2,nsamp,ampst=0.8,amped=1.2){
    res= cmdscale(dist2, k = 2,eig=TRUE);
    eig=res$eig;
    p1 = sum(abs(eig[1:2])) / sum(abs(eig)) 
    p2 = sum(eig[1:2] ^ 2) / sum(eig ^ 2) 
    cat("Mardia fit measure 1 = ", p1, "\n")
    cat("Mardia fit measure 2 = ", p2, "\n")
    
    dis2=(dist(res $points));
    dis2min=min(dis2);
    dis2max=max(dis2);
    dis2=as.matrix(dis2);
    
    edge_st=dis2min*ampst;
    edge_ed=dis2max*amped;
    edges=exp(seq(log(edge_st),log(edge_ed),,nsamp));

    
    return(list(points=res$points,dis2=dis2,p1=p1,p2=p2,edges=edges));
}

##########################################################
cerwave <- function(i,edge){
    if((gidbuf[i]==0)){
        gidbuf[i] <<- ngroupbuf;
        lis=which(distbuf[i,] < edge);
 
        if(length(lis)>0){
            for(j in 1:length(lis))cerwave(lis[j],edge);
        }
    }
}


gidbuf <<- 0;
ngroupbuf <<- 1;
distbuf <<- 0.0;
do_cerwave <- function(dis2,edge){
    distbuf <<- dis2;
    gidbuf <<- rep(0,nrow(distbuf));
    ngroupbuf <<- 1;      
    for(i in 1:length(gidbuf)){
        if(gidbuf[i]==0){     
            cerwave(i,edge);
            ngroupbuf <<- ngroupbuf+1;
        }        
    }    
    ngroupbuf <<- ngroupbuf-1;    
    return(list(ngroup=ngroupbuf,gid=gidbuf));
}


calc_groups <- function(dist2,edges,param=param0){
    if(param$flag_fix_kdis==1){
        kdis=param$kdis0;
    }else{
        res= cmdscale(dist2, k = 2,eig=TRUE);
        eig = res$eig                        
        cum_eig=cumsum(eig^2)/sum(eig^2);
        kdis=max(2,min(which(cum_eig>0.95)));
    }
    
    print("kdis:");
    print(kdis);
    resd= cmdscale(dist2, k = kdis,eig=TRUE);
    dis2=(dist(resd$points));
    dis2min=min(dis2);
    dis2max=max(dis2);
    dis2=as.matrix(dis2);
    
    dismin=min(dist2[which(dist2 > 0)]);
    
    
    print("cerwave..");
    nn=length(edges);
    ##ngroups2=rep(0,length(edges));    
    ##for(i in 1:length(edges)){
    ngroups2=foreach::foreach(i=1:nn,.combine='rbind',.packages="foreach")%do%{
        if(i%%50==0)cat(100*i/nn,"% ");
                
        edge = edges[i];
        re_cerwave=do_cerwave(dis2,edge);
        ngroup=re_cerwave$ngroup;
        gid=re_cerwave$gid;
        ngroup;
    };
    
    print("done");
    return(ngroups2);
}

calc_groups_sm <- function(dist2,edges,edge2ampsm,ndiv=100,flag_pforeach=param0$flag_pforeach){
    res= cmdscale(dist2, k = 2,eig=TRUE);
    ##ngroups_sm=rep(0,length(edges));    
    ##for(i in 1:length(edges)){
    ##    ngroups_sm[i]=(find_peak(res,edges[i]*edge2ampsm,ndiv=ndiv))$ngroup;
    ##}
    nn=length(edges);
 

    if(flag_pforeach==1){
        print("using pforeach..");
        ngroups_sm=pforeach::pforeach(i=1:nn, .combine='rbind',.errorhandling="stop")({    
            if(i%%50==0)cat(100*i/nn,"% ");
            (find_peak(res,edges[i]*edge2ampsm,ndiv=ndiv))$ngroup; 
        });
    }else{
        print("using foreach..");
        ngroups_sm=foreach::foreach(i=1:nn,.combine='rbind',.packages="foreach")%do%{
            if(i%%50==0)cat(100*i/nn,"% ");
            (find_peak(res,edges[i]*edge2ampsm,ndiv=ndiv))$ngroup;
            
        };
    }

    
    print("done");
    
    ngroups_sm_ob=sort(unique(ngroups_sm));
    ampsm_ob=ngroups_sm_ob*0.0;
    for(i in 1:length(ngroups_sm_ob)){
        ampsm_ob[i]=edge2ampsm*mean(edges[ngroups_sm==ngroups_sm_ob[i]]);
    }
    return(list(ngroups_sm=ngroups_sm,ngroups_sm_ob=ngroups_sm_ob,ampsm_ob=ampsm_ob));
}


czz <<- 0.0;
cndiv <<- 0;
ccerwave <- function(i,j,id){
    ccer[i,j,id] <<- 1;
    if((i<cndiv)&&(ccer[i+1,j,id]==0)&&(czz[i+1,j]<czz[i,j]))ccerwave(i+1,j,id);
    if((i>1)&&(ccer[i-1,j,id]==0)&&(czz[i-1,j]<czz[i,j]))ccerwave(i-1,j,id);
    if((j<cndiv)&&(ccer[i,j+1,id]==0)&&(czz[i,j+1]<czz[i,j]))ccerwave(i,j+1,id);
    if((j>1)&&(ccer[i,j-1,id]==0)&&(czz[i,j-1]<czz[i,j]))ccerwave(i,j-1,id);

    if((i<cndiv)&&(j<cndiv)&&(ccer[i+1,j+1,id]==0)&&(czz[i+1,j+1]<czz[i,j]))ccerwave(i+1,j+1,id);
    if((i<cndiv)&&(j>1)&&(ccer[i+1,j-1,id]==0)&&(czz[i+1,j-1]<czz[i,j]))ccerwave(i+1,j-1,id);
    if((i>1)&&(j<cndiv)&&(ccer[i-1,j+1,id]==0)&&(czz[i-1,j+1]<czz[i,j]))ccerwave(i-1,j+1,id);
    if((i>1)&&(j>1)&&(ccer[i-1,j-1,id]==0)&&(czz[i-1,j-1]<czz[i,j]))ccerwave(i-1,j-1,id);
    
}

ccerwaved <- function(i,j,id,wa){
    ccer[i,j,id] <<- wa;
    wd=1;
    if((i<cndiv)&&(ccer_fix[i+1,j,id]==0)&&(ccer[i+1,j,id]>wa+wd)&&(czz[i+1,j]<czz[i,j]))ccerwaved(i+1,j,id,wa+wd);
    if((i>1)&&(ccer_fix[i-1,j,id]==0)&&(ccer[i-1,j,id]>wa+wd)&&(czz[i-1,j]<czz[i,j]))ccerwaved(i-1,j,id,wa+wd);
    if((j<cndiv)&&(ccer_fix[i,j+1,id]==0)&&(ccer[i,j+1,id]>wa+wd)&&(czz[i,j+1]<czz[i,j]))ccerwaved(i,j+1,id,wa+wd);
    if((j>1)&&(ccer_fix[i,j-1,id]==0)&&(ccer[i,j-1,id]>wa+wd)&&(czz[i,j-1]<czz[i,j]))ccerwaved(i,j-1,id,wa+wd);

    wd=sqrt(2);
  if((i<cndiv)&&(j<cndiv)&&(ccer_fix[i+1,j+1,id]==0)&&(ccer[i+1,j+1,id]>wa+wd)&&(czz[i+1,j+1]<czz[i,j]))ccerwaved(i+1,j+1,id,wa+wd);
    if((i<cndiv)&&(j>1)&&(ccer_fix[i+1,j-1,id]==0)&&(ccer[i+1,j-1,id]>wa+wd)&&(czz[i+1,j-1]<czz[i,j]))ccerwaved(i+1,j-1,id,wa+wd);
    if((i>1)&&(j<cndiv)&&(ccer_fix[i-1,j+1,id]==0)&&(ccer[i-1,j+1,id]>wa+wd)&&(czz[i-1,j+1]<czz[i,j]))ccerwaved(i-1,j+1,id,wa+wd);
    if((i>1)&&(j>1)&&(ccer_fix[i-1,j-1,id]==0)&&(ccer[i-1,j-1,id]>wa+wd)&&(czz[i-1,j-1]<czz[i,j]))ccerwaved(i-1,j-1,id,wa+wd);
    
##    if((i<cndiv)&&(j<cndiv)&&(ccer_fix[i+1,j+1,id]==0)&&(ccer[i+1,j+1,id]>wa+1)&&(czz[i+1,j+1]<czz[i,j]))ccerwaved(i+1,j+1,id,wa+wd);
##    if((i<cndiv)&&(j>1)&&(ccer_fix[i+1,j-1,id]==0)&&(ccer[i+1,j-1,id]>wa+1)&&(czz[i+1,j-1]<czz[i,j]))ccerwaved(i+1,j-1,id,wa+wd);
##    if((i>1)&&(j<cndiv)&&(ccer_fix[i-1,j+1,id]==0)&&(ccer[i-1,j+1,id]>wa+1)&&(czz[i-1,j+1]<czz[i,j]))ccerwaved(i-1,j+1,id,wa+wd);
##    if((i>1)&&(j>1)&&(ccer_fix[i-1,j-1,id]==0)&&(ccer[i-1,j-1,id]>wa+1)&&(czz[i-1,j-1]<czz[i,j]))ccerwaved(i-1,j-1,id,wa+wd);
    
}

map_power <- function(pcoa,rec,ampsm_opt,ndiv=100){
    if(flag_envs).ee.append("map_power",environment());

    popids2=rec$popids2;
   
    px = pcoa$points[,1]
    py = pcoa$points[,2]
   
    ampsm = ampsm_opt;
    
    res_peak=find_peak(pcoa,ampsm_opt,ndiv=ndiv);
    peak=res_peak$peak;
    lispeak=res_peak$lisp;
    peakx=res_peak$peakx;
    peaky=res_peak$peaky;
    xx0=res_peak$xx0;
    yy0=res_peak$yy0;
    xx=res_peak$xx;
    yy=res_peak$yy;
    sx=res_peak$sx;
    sy=res_peak$sy;
    zz = res_peak$zz;
    czz <<- res_peak$zz;
    cndiv <<- ndiv;
    
    print(peak);

    ngroup = length(lispeak);
    print(ngroup);

    ccer <<- array(rep(as.integer(px)*0,ngroup),dim=c(ndiv,ndiv,ngroup));

    for(k in 1:(length(lispeak))){
             ccerwave(peak[1,k],peak[2,k],k);
    }
    
    for(i in 1:ndiv){
        for(j in 1:ndiv){
            if(sum((ccer[i,j,]>0))>1){
                ccer[i,j,] <<- 0;
            }        
        }
    }
    
    ccer0 = ccer;

    ccer_fix <<- ccer;
    ccer[which(ccer==0)] <<- 100000;
    for(k in 1:(length(lispeak))){
        for(i in 1:ndiv){
            for(j in 1:ndiv){
                if(ccer_fix[i,j,k]==1)ccerwaved(i,j,k,1);
            }
        }
    }
    
    ccer_buf = ccer;
    ccer_buf[which(ccer_buf==0)]=100000;
    ccer1 = ccer*0;
    for(i in 1:ndiv){
        for(j in 1:ndiv){
            mid=which.min(ccer_buf[i,j,]);
            ccer1[i,j,] = 0;
            ccer1[i,j,mid] = 1;        
        }
    }
    
    cmask = array(rep(0,length(px)*ngroup),dim=c(length(px),ngroup));
    for(k in 1:ngroup){
        cmask[,k] = interp2(yy0,xx0,ccer1[,,k],py,px);
    }
    cid_sm = rep(0,length(px));
    for(i in 1:length(px)){
        cid_sm[i]=which.max(cmask[i,]);
    }

##    print("AAA");

  ##  print(lispeak);
    return(list(cid_sm=cid_sm,peakx=peakx,peaky=peaky,xx0=xx0,yy0=yy0,xx=xx,yy=yy,zz=zz,sx=sx,sy=sy,ccer0=ccer0,ccer1=ccer1,ndiv=ndiv));
}


pfunc <- function(px,py,xx0,yy0,peakx,peaky,ccer0,ccer1,cid_sm,labels,lev=0.95,pals= NULL,flag_map=1){
    ngroup=max(cid_sm);
    if(length(pals)==0)pals=rainbow(ngroup,v=1.0);
    
    if(flag_map==0){
        text(x=px,y=py,col=pals[cid_sm],labels=labels,cex=0.7,font=1)
    }else{
        
        
        for(ii in 1:ngroup){
            ##        points(xx[which(ccer1[,,ii]>0)],yy[which(ccer1[,,ii]>0)],col=pals[ii],pch=20,cex=0.7);
            contour(xx0,yy0,ccer1[,,ii],add=TRUE,drawlabels=FALSE,method="simple",col=pals[ii],levels=c(lev),labels=NULL,lwd=2,lty=1)
            contour(xx0,yy0,ccer0[,,ii],add=TRUE,drawlabels=FALSE,method="simple",col=pals[ii],levels=c(lev),labels=NULL,lwd=1,lty=2)   
            ##        text(xx[which(ccer1[,,ii]>0)],yy[which(ccer1[,,ii]>0)],label=paste(ii),col=pals[ii],pch=20,cex=0.5);
        }
        
        text(x=px,y=py,col=pals[cid_sm],labels=labels,cex=0.7,font=1)
        
        points(peakx,peaky,col="black",pch=17,lwd=3,cex=1.8);
        points(peakx,peaky,col=pals,pch=17,lwd=3,cex=1.5);
        points(peakx,peaky,col="white",pch=17,lwd=3,cex=0.8);
    }
    
}


sqm <-function(x){
    return(sum((x-mean(x))^2));
}

calc_inde <- function(pxp,pyp){
    vara=sqm(pxp)+sqm(pyp);
    varw=0.0;
    varb=0.0;
    for(k in 1:ngroup){
        px1=pxp[which(cid_sm==k)];
        py1=pyp[which(cid_sm==k)];
        varw= varw+ sqm(px1)+sqm(py1);
        varb=varb+length(px1)*((mean(px1)-mean(pxp))^2+(mean(py1)-mean(pyp))^2);
       ## varw= varw+ sqm(px1);
       ## varb=varb+length(px1)*((mean(px1)-mean(pxp))^2);
    }

##    vb=varb/(ngroup-1);
    vw=varw/(length(pxp)-ngroup);
    vb=varb/(length(pxp)-ngroup);
    inde=(vb/vw);
    ##inde=(varb/varw);
    ##inde=(varb/vara);
    return(list(inde=inde,vb=(varb/(length(pxp))),va=vara,vw=(varw/(length(pxp)-ngroup))));
}

cmdist <- function(dist2p,cid_sm,test=FALSE,flag_welch_t_test=param0$flag_welch_t_test,nperm=param0$nperm,flag_pforeach=param0$flag_pforeach){
    
    ngroup=max(cid_sm);
    gdis=matrix(rep(0.0,ngroup*ngroup),nrow=ngroup,ncol=ngroup);
    pvs=gdis;

    ilis=NULL;
    jlis=NULL;
    for(i in 1:(ngroup-1)){
        for(j in (i+1):ngroup){
            ilis=c(ilis,i);
            jlis=c(jlis,j);
        }
    }


if(0){
    for(k in 1:length(ilis)){
        res=cmdist_unit(k, ilis,jlis,dist2p,cid_sm,test=test,flag_welch_t_test=flag_welch_t_test,nperm=nperm);

    }
}

     if(flag_pforeach==1){
         print("using pforeach..");
         res=pforeach::pforeach(k=1:length(ilis), .combine='rbind',.errorhandling="stop")({    
             cmdist_unit(k, ilis,jlis,dist2p,cid_sm,test=test,flag_welch_t_test=flag_welch_t_test,nperm=nperm);
         });
         
     }else{
         print("using foreach..");
         res=foreach::foreach(k=1:length(ilis),.combine='rbind',.packages="foreach")%do%{
             cmdist_unit(k, ilis,jlis,dist2p,cid_sm,test=test,flag_welch_t_test=flag_welch_t_test,nperm=nperm);
         };

     }
         
##    print(c(length(ilis),res));

    for(k in 1:length(ilis)){
        if(length(ilis)==1){
            dm0=res[1];
            pv0=res[2];
            
        }else{
            dm0=res[k,1];
            pv0=res[k,2];
            
        }
        
        gdis[ilis[k],jlis[k]]=dm0;
        gdis[jlis[k],ilis[k]]=dm0;
        pvs[ilis[k],jlis[k]]=pv0;
        pvs[jlis[k],ilis[k]]=pv0;
    }
    ##}
    

    ##    print("a");
    ##    print(pvs);
    return(list(gdis=gdis,pvs=pvs));
}

#' @export
cmdist_unit <- function(k, ilis,jlis,dist2p,cid_sm,test=FALSE,flag_welch_t_test=1,nperm=1000){
    i=ilis[k];
    j=jlis[k];
    lis=which((cid_sm==i)+(cid_sm==j)>0);
    dist2pp=dist2p[lis,lis];
    cid=cid_sm[lis];
    rec=cmdscale(dist2pp);
    pp=rec[,1];
    pi=(rec[(cid==i),1]);
    pj=(rec[(cid==j),1]);
    sdi=sd(pi);
    sdj=sd(pj);
    dij=abs(mean(pi)-mean(pj));
    pv=1.0;
    ##sdw=(length(pi)*sdi+length(pj)*sdj)/(length(lis));
    sdw=sqrt((sqm(pi)+sqm(pj))/(length(lis)-2));
    if(flag_welch_t_test==1)sdw1=sqrt(sdi^2/length(pi)+sdj^2/length(pj));
    
    ##gdis[i,j]=dij/(0.5*(sdi+sdj));
    ##gdis[i,j]=dij/sdw;
    ##gdis[j,i]=gdis[i,j];
            dm0=dij/sdw;
            
            if(test==TRUE){
                if(flag_welch_t_test==1)sdw=sdw1;
                inde0=dij/sdw;
                
                
                ##nperm=1000;
                indes=rep(0.0,nperm);
                for(k in 1:nperm){
                    lis=order(rnorm(length(pp)));
                    pi=(rec[(cid[lis]==i),1]);
                    pj=(rec[(cid[lis]==j),1]);
                    sdi=sd(pi);
                    sdj=sd(pj);
                    dij=abs(mean(pi)-mean(pj));
                    sdw=sqrt((sqm(pi)+sqm(pj))/(length(lis)-2));
                    if(flag_welch_t_test==1)sdw1=sqrt(sdi^2/length(pi)+sdj^2/length(pj));
                    if(flag_welch_t_test==1)sdw=sdw1;
                    indes[k]=dij/sdw;
                }
                p_value=(sum(indes>inde0)+sum(indes==inde0)*0.5)/nperm;

##                print(c(dm0,p_value));
                
                if(!is.na(p_value)){
                     pv=p_value;
                     ##pvs[j,i]=p_value;
                }
            }
   
            return(c(dm0,pv))
     }


##' @title popmap.plot". 
#' @export
popmap.plot <- function(popmap,param=param0,perm=0,resp=NULL,main=NULL){
    if(flag_envs).ee.append("popmap.plot",environment());
    cid_sm=popmap$cid_sm;
    if(perm==1){
        pcoap=popmap$pcoap;
        powermap=popmap$powermap;

        if(length(resp)==0)resp=find_peak(pcoap,popmap$ampsm,ndiv=powermap$ndiv);

        if(length(main)==0){
            if(popmap$maxp>0){
                maxij=(which(popmap$gdist$pvs==max(popmap$gdist$pvs),arr.ind=T))[1,];
                main=sprintf("max p-value: %f (%d:%d nperm: %d)",max(popmap$gdist$pvs),maxij[1],maxij[2],param$nperm);
            }else{
                main=sprintf("max p-value: %f (nperm: %d)",popmap$maxp,param$nperm);
             
            }
            
        }



        plot_power_map(pcoap,resp,cid_sm,flag_fix_lim=param$flag_fix_lim,pals=popmap$pals,ccer0=powermap$cerr0,ccer1=powermap$ccer1,labels=paste(cid_sm),main=main,flag_map=0,param=param);

        
    }else{
        pcoa=popmap$pcoa;
        powermap=popmap$powermap;
        if(length(main)==0)main=sprintf("ngroups: %d   sd: %2.3f\n  max p-value: %f    nperm: %d",popmap$ngroup, popmap$ampsm,popmap$maxp,param$nperm);

        plot_power_map(pcoa,powermap,cid_sm,flag_fix_lim=0,pals=popmap$pals,ccer0=powermap$ccer0,ccer1=powermap$ccer1,labels=popmap$popids2,flag_map=1,main=main,param=param);
        
        
    }
}

plot_power_map <- function(pcoa,res,cid_sm,flag_fix_lim=0,pals=NULL,ccer0=NULL,ccer1=NULL,labels=NULL,flag_map=1,main=NULL,param=param0){
    if(flag_envs).ee.append("plot_power_map",environment());
    xx0=res$xx0;
    yy0=res$yy0;
    xx=res$xx;
    yy=res$yy;
    zz=res$zz;
    sx=res$sx;
    sy=res$sy;

    
    ngroup = max(cid_sm);
    peakx = res$peakx;
    peaky = res$peaky;
    
    px=pcoa$points[,1];
    py=pcoa$points[,2];
    if(length(pals)==0)pals = rainbow(ngroup,v=1.0);

        
    par(mar=c(5,5,5,6));
    par(xpd=T);
    dev.hold();
    
    if(length(pals)==0)pals = rainbow(ngroup,v=1.0);
    
    
    nlev=20;
    
    if(param$flag_fix_lim==0){
        plot(xx0,yy0,type="n",xlim=c(min(xx0),max(xx0)),ylim=c(min(yy0),max(yy0)),main=main,xaxs="i",yaxs="i",xlab="Axis 1", ylab="Axis 2");
    }else{    
        plot(xx0,yy0,type="n",xlim=xlim0,ylim=ylim0,main=main,xaxs="i",yaxs="i",xlab="Axis 1", ylab="Axis 2");
        
        polygon(c(xlim0[1],xlim0[2],xlim0[2],xlim0[1]),c(ylim0[1],ylim0[1],ylim0[2],ylim0[2]), col=(gray.colors(nlev,start=0.1,end=0.6))[1]);
        
    }
    

    .filled.contour(xx0,yy0,zz,levels=seq(0,max(zz),,nlev),col=gray.colors(nlev,start=0.1,end=0.6));

    
    pfunc(px,py,xx0,yy0,peakx,peaky,ccer0,ccer1,cid_sm,labels,pals=pals,flag_map=flag_map);
    
    legend(par()$usr[2], par()$usr[4], legend = paste("group",seq(length(pals))), col = pals,pch=rep(15,length(pals)),cex=0.8,text.col="white",bg="#555555");

        dev.flush();
}


print_gdist <- function(gdist,gdist0){
##                options(digits=4);
                
                cat("\n\nGenetic distance matrix calculated from data for grouping: \n");
                print(gdist0$gdis);
                
                cat("\n\nGenetic distance matrix calculated from data for petmutation test: \n");
                print(gdist$gdis);
##                options(digits=7);     
                cat("p-values: \n");
                print(gdist$pvs);
                cat("\n");
}

plot_power <- function(pcoa,rec,ampsm_opt,param=param0){
    if(flag_envs).ee.append("plot_power",environment());
    dist2=rec$dist2;
    dist2p=rec$dist2p;
    popids2=rec$popids2;

    ndiv=param$ndiv;
    
    res=map_power(pcoa,rec,ampsm_opt,ndiv=param$ndiv);
    cid_sm = res$cid_sm;
    ngroup = max(cid_sm);


    dev.set(4);

    main=sprintf("sd: %2.3f     groups: %d",ampsm_opt,ngroup);

    pals = rainbow(ngroup,v=1.0);


        
        plot_power_map(pcoa,res,cid_sm,flag_fix_lim=param$flag_fix_lim,pals=pals,ccer0=res$ccer0,ccer1=res$ccer1,labels=popids2,flag_map=1,main=main,param=param);


        
    dev.set(5);
    
    
    if(abs(param$flag_perm)==1){
        cid_count = rep(0,max(cid_sm));
        for(i in 1:max(cid_sm)){
            cid_count[i] = sum(cid_sm==i);
        }
        
        flag_do_perm=0;
        if(sum(cid_count<=1)<2)flag_do_perm=1;

        gdist=NULL;
        gdist0=NULL;
        if(flag_do_perm==1){
            print("calculating genetic distances among subpopulations..");
            gdist=cmdist(dist2p,cid_sm,test=TRUE,flag_welch_t_test=param$flag_welch_t_test,nperm=param$nperm,flag_pforeach=param$flag_pforeach);
            
            gdist0=cmdist(dist2,cid_sm,test=FALSE,flag_welch_t_test=param$flag_welch_t_test,nperm=param$nperm,flag_pforeach=param$flag_pforeach);
        }

        
        print("done");

        pcoap=cmdscale(dist2p,k=2,eig=TRUE);

        px=pcoa$points[,1];
        py=pcoa$points[,2];
        pxp=pcoap$points[,1];
        pyp=pcoap$points[,2];
        
        pxp=sign(sum(pxp*px))*pxp;
        pyp=sign(sum(pyp*py))*pyp;
        pcoap$points[,1]=pxp;
        pcoap$points[,2]=pyp;
        
        
        resp=find_peak(pcoap,ampsm_opt,ndiv=ndiv);
        

        maxij=1.0;
        if(flag_do_perm==1){
            maxij=(which(gdist$pvs==max(gdist$pvs),arr.ind=T))[1,];
        }    
        if(flag_do_perm==1){
            main0=sprintf("max p-value: %f (%d:%d nperm: %d)",max(gdist$pvs),maxij[1],maxij[2],param$nperm);
        }else{
            main0="too small populations!!";
        }

        
        plot_power_map(pcoap,resp,cid_sm,flag_fix_lim=param$flag_fix_lim,pals=pals,ccer0=res$cerr0,ccer1=res$ccer1,labels=paste(cid_sm),main=main0,flag_map=0,param=param);

        if(flag_do_perm==1){

            print_gdist(gdist,gdist0);
        }
        
        if(0){
            for(i in 1:ngroup){
                cat(i,":\n ");
                for(j in 1:ngroup){
                    if(i!=j){
                        cat(sprintf("%d %2.2f ",j,gdis[i,j]));
                    }
                }
                cat("\n");
            }
        }

        ##hist(indes,xlim=c(0,inde0*1.2));
        ##points(inde0,0,pch=17,col="red");
}
    maxp=1.0;
      if(flag_do_perm==1){
          maxp=max(gdist$pvs);
      }
    
    return(list(ngroup=ngroup,maxp=maxp,cid_sm=cid_sm,popids2=popids2,gdist=gdist,gdist0=gdist0,powermap=res,pcoa=pcoa,pcoap=pcoap,ampsm=ampsm_opt,pals=pals));
}


pngout<-function(dev_id=as.numeric(dev.list()[length(dev.list())]),plotfile="testout",density=150,geometry=600,outeps=FALSE,prefix="./",show=TRUE,outpng=TRUE){
    
    dens=as.character(density);
    geom=as.character(as.integer(geometry));
    dev.set(dev_id);

    tmpf=paste0("R_pngout_temp");
    
    dev.copy2eps(file=sprintf("%s.eps",tmpf));
   
     if(outeps){
         ##system(paste("cp .R_pngout_temp.eps ",prefix,plotfile,".eps",sep=""));
         system(sprintf("cp %s.eps %s%s.eps",tmpf,prefix,plotfile));
         cat("eps output:",paste(prefix,plotfile,".eps",sep=""),"\n")
     }
    
    if(outpng){
        system(sprintf("convert -density %sx%s -geometry %s  -background white -alpha remove %s.eps %s.png",dens,dens,geom,tmpf,tmpf));
        system(sprintf("mv %s.png %s%s.png",tmpf,prefix,plotfile));
        cat(sprintf("png output: %s%s.png \n",prefix,plotfile));
        if(show)system(sprintf("display %s%s.png&",prefix,plotfile));
    
    }
    system(sprintf("rm %s.eps",tmpf));
    
}

##' @title popmap.plotdata". 
#' @export
popmap.plotdata <- function(pcoa,data,param=param0,sample_group_name=NULL,flag_lim_out=0){
   if(flag_envs).ee.append("popmap.plotdata",environment());

   if(length(sample_group_name)==0)sample_group_name=data$sample_group_name;
   
   amp_euc=data$amp_euc;
   amp_bin=data$amp_bin;
   omit_popid=data$omit_popid;
   
   runnames=data$runnames;
    masks=data$masks;
    mask=data$mask;
    nrun=data$nrun;
    popnames=data$popnames;
    popids2=data$popids2;
    popids=data$popids;
    runids1=data$runids1;    
   runid=data$runid;
   
    n_eu=data$neu;
    nn=nrow(data$tab1);
    nc=ncol(data$tab1);
    xmax=max(pcoa$points[,1]);
    xmin=min(pcoa$points[,1]);
    gx=0.5*(xmax+xmin);
    dx=0.5*(xmax-xmin);
    ymax=max(pcoa$points[,2]);
    ymin=min(pcoa$points[,2]);
    gy=0.5*(ymax+ymin);
    dy=0.5*(ymax-ymin);

    p1=pcoa$p1;
    p2=pcoa$p2;

   rmask=(rowSums(masks)>0);

   
xlim0=1.2*dx*c(-1,1);
ylim0=1.2*dy*c(-1,1);
xlim0=xlim0+gx;
ylim0=ylim0+gy;

    
    ##    tmar=sprintf("\n Edge_dis: %2.4f, dimdis: %d, Mardia-fit measure 1: %2.2f, 2:%2.2f",edge,kdis,p1,p2);
        tmar=sprintf("\n Mardia-fit measure 1: %2.2f, 2:%2.2f",p1,p2); 

    if(param$edge_null>1){
        edge_nullt="None";
    }else{
        edge_nullt=param$edge_null;
    }
    
    if(param$flag_binary==0)main_name=paste(sample_group_name,"(euclid-distance)",", null edge:",edge_nullt,", loci:",nc,tmar);
if(param$flag_binary==1)main_name=paste(sample_group_name,"(my-binary-distance)",", null edge:",edge_nullt,", loci:",nc,tmar);

   ##if(param$flag_binary==2)main_name=paste(sample_group_name,"(hybrid-distance)",", null edge:",edge_nullt,", loci:",nc,"\n amp_euc:",amp_euc,", amp_bin:",amp_bin,", mean non-NULL pairs:", as.integer(mean(n_eu)),tmar);
   if(param$flag_binary==2)main_name=paste(sample_group_name," , loci:",nc,"\n amp_euc:",amp_euc,", amp_bin:",amp_bin,", mean non-NULL pairs:", as.integer(mean(n_eu)),tmar);
if(param$flag_binary==3)main_name=paste(sample_group_name,"(R-binary-distance)",", null edge:",edge_nullt,", loci:",nc,tmar);

 ## 0 euclid, 1 my binary, 2 hybrid, 3 binary by R function

par(mar=c(5,5,5,6));
par(xpd=T);
if(param$flag_runid==1){
    plot.new();    
    ##plot(pcoa,type="n",main=paste0(main_name,"\n run_id: color,   population_id: number"));
    plot(pcoa$points,type="n",main=paste0(main_name),xlab="Axis 1",ylab="Axis 2",xlim=xlim0,ylim=ylim0,cex.main=0.8);
##m12: 3,10,21,20,9,5,11,1,13,8
##m13: 22,13,7,12,19,4,6,5,9,14,10,
##m14: 1,3,4,5,6,7,8,9,10,11,12,13,14,18,20,21,22

##pal=c("yellow","orange","green","blue","purple","red");
##pal=rev(rainbow(6,end=0.7));

    if(param$flag_cols_cerwave==1){
        ##        pal=rev(rainbow(ngroup,end=0.7));
        pal=rev(rainbow(ngroup,end=0.7,v=0.8));
        cols=pal[gid];
        
    }else{
        ##        pal=rev(rainbow(nrun,end=0.7));
        ##        pal=rev(rainbow(nrun,end=0.7,v=0.8));
        runid=data$runid;
        pal=rev(rainbow(length(runid),end=0.7,v=0.8));
        cols=rep("blue",nn);
        
        ##for(i in 1:nrun)cols[masks[i,]]=pal[i];
        for(i in 1:length(runid))cols[masks[runid[i],]]=pal[i];
    }
    
    
    ##text(x=pcoa$points[mask,1],y=pcoa$points[mask,2],col=cols[mask],labels=popids2,cex=0.7,font=2)

    
    text(x=pcoa$points[mask,1],y=pcoa$points[mask,2],col=cols[mask],labels=popids2,cex=0.7,font=1)

    run_name="run";
    if(param$flag_cols_cerwave==1)run_name="group";


    
    ##legend(param$legend_side, legend = paste(run_name,seq(length(pal))[rmask]), col = pal[rmask],pch=rep(15,length(pal[rmask])),cex=0.8);

    legend(param$legend_side, legend = paste(run_name,data$runnames), col = pal,pch=rep(15,length(pal)),cex=0.8);
    
##    legend(par()$usr[2], par()$usr[3], legend = paste("run",seq(length(pal))), col = pal,pch=rep(15,length(pal)))

if(param$flag_show_legend_location==1){
if(param$flag_read_sample_location==1){
    ##popid_u=unique(popids);
    ##popname_u=rep("",length(popid_u));
    ##for(i in 1:length(popid_u))popname_u[i]=popnames[(which(popids==popid_u[i]))[1]];
    
    
    legend(par()$usr[2], par()$usr[4], legend = paste0(data$popid.legend,": ",data$popname.legend), col = "black",cex=0.7,bty="n");

    text(par()$usr[2], par()$usr[3], labels = paste("Omitted \npopids:",paste(omit_popid,collapse=", ")), col = "red",cex=0.8,pos=4)

    ##mtext(paste0("Omitted popids:",omit_popid), side=4,col = "red",cex=0.8)
    }
}
}

if(param$flag_runid==0){
    plot.new();    
    plot(pcoa$points,type="n",main="run_id: number,   population_id: color",xlab="Axis 1");

    pops=sort(unique(popids));
    pal=rev(rainbow(length(pops),end=0.7));
    pal=rev(rainbow(length(pops),end=0.7,v=0.8));
    
    cols=rep("blue",nn);
    for(i in 1:length(pops)){
        lis=which(popids2==pops[i]);
        cols[lis]=pal[i];
    }
##    mask=mask1;
    

    text(x=pcoa$points[mask,1],y=pcoa$points[mask,2],labels=runids1[mask],cex=0.7)
    text(x=pcoa$points[mask,1],y=pcoa$points[mask,2],col=cols[mask],labels=runids1[mask],cex=0.7)

    legend("topright", legend = pops, col = pal,pch=rep(3,length(pal)));
}

   if(flag_lim_out==1)return (list(xlim0=xlim0,ylim0=ylim0));
}


##' @title popmap.find". 
#' @export
popmap.find <- function(rec,param=param0,sample_group_name=NULL){
    if(flag_envs).ee.append("popmap.find",environment());

    if(length(sample_group_name)==0)sample_group_name=rec$sample_group_name;
    if(length(dev.list())==0)xinit(param$win_width,param$win_height);
    
    dist2=rec$dist2;
    dist2p=rec$dist2p;    
    pcoa=get_pcoa2d(dist2,param$nsamp);
    edges=pcoa$edges;
    
    print("find groups for different edge values..");
    ngroups2=calc_groups(dist2,edges,param=param);
    
    print("find groups (peaks) for different gaussian smoothing..");
    re_calc_groups_sm=calc_groups_sm(dist2,edges,param$edge2ampsm,ndiv=param$ndiv,flag_pforeach=param$flag_pforeach);
    ngroups_sm=re_calc_groups_sm$ngroups_sm;
    ampsm_ob=re_calc_groups_sm$ampsm_ob;
    ngroups_sm_ob=re_calc_groups_sm$ngroups_sm_ob;
    
############################################################################
    dev.set(2);

lims=popmap.plotdata(pcoa,rec,param=param,sample_group_name=sample_group_name,flag_lim_out=1);
xlim0=lims$xlim0;
ylim0=lims$ylim0;
#####################################################################

dev.set(3);
    edge2ampsm=param$edge2ampsm;
    edgebuf=c(edges,edges*edge2ampsm);
plot(edges*edge2ampsm,ngroups_sm,type="l",log="xy",col="blue",ylim=c(1,ncol(dist2)),xlim=c(min(edgebuf),max(edgebuf)),xlab="sd for Gaussian smoothing",ylab="Number of groups")
lines(edges,ngroups2,col="black");

    legend("topright", legend = c("Gaussian smooth","Neighbor connect","Significance"), col = c("blue","black","red"),pch=c(15,15,1),cex=0.8);

    popms=list(pcoa=pcoa);
    if(param$ngroup_show<0){
        pvs=rep(1.0,20);
        count_sig=0;
        for(k in 2:20){
            cat("\n\n::::Testing significance for ",k, "groups-mapping  (ampsm:",ampsm_ob[k],") \n");
            maxp=1.0;
            popm=plot_power(pcoa,rec,ampsm_ob[k],param=param);
            maxpb=popm$maxp;
            
            if(!is.na(maxpb))maxp=maxpb;
            pvs[k]=maxp;
            
            if(maxp<param$edge_p_value){
            count_sig=count_sig+1;
            ##            si=data.frame(rec$sample.info,GroupID=rep(0,nrow(rec$sample.info)),pcoa.x=pcoa$points[,1],pcoa.y=pcoa$points[,2]);
            buf=rep(0,nrow(rec$sample.info));
            si=data.frame(rec$sample.info,GroupID=buf,pcoa.x=buf,pcoa.y=buf);
            sname=rec$sample.name2;
            ngroup=popm$ngroup;
            for(i in 1:nrow(si)){
                id=which(sname==si$SampleName2[i]);
                if(length(id)>0){
                    si$GroupID[i]=popm$cid_sm[id];
                    si$pcoa.x[i]=pcoa$points[id,1];
                    si$pcoa.y[i]=pcoa$points[id,2];
                    
                }
            }
##            sample.info=si;
            nf=which(colnames(si)=="SampleName");
            ##sample.info=si[,c(seq(2,nf),1,seq(nf+1,ncol(si)))];
            sample.info=si[,c(seq(2,ncol(si)),1)];
            popm=append(popm,list(sample.info=sample.info));

            popms=append(list(popm),popms);
            names(popms)[1]=paste0("ngroup",ngroup);
            
            dev.set(3);
            
            points(ampsm_ob[k],ngroups_sm_ob[k],pch=1,col="red",lwd=1);
            dev.set(4);
        }else{
            print(sprintf("p-value %2.3f higher than %2.3f  ->>  not sifnificant !!!",maxp,param$edge_p_value));
            
            if((ngroups_sm_ob[k]>param$group_max)||(k==length(ngroups_sm_ob))){
                if(count_sig>0){
                    ##opid=max(which(pvs<param$edge_p_value));
                    ##popm=plot_power(pcoa,rec,ampsm_ob[opid],param=param);
                    
                    dev.set(4);
                    popmap.plot(popms[[1]],param=param,perm=0)
                    dev.set(5);
                    popmap.plot(popms[[1]],param=param,perm=1)
                    
                    ##maxpb=popm$maxpb;
                    cat("\n\n::::Maximum identified subpopulations: ",popms[[1]]$ngroup,"\n");
                    cat("Max p-value: ",popms[[1]]$maxp,"\n");

                    print_gdist(popms[[1]]$gdist,popms[[1]]$gdist0);



                ##popms=append(list(popm),popms);
                ##names(popms)[1]=paste0("ngroup",ngroup);
                print(names(popms));
                }else{
                    print("No significant group structure detected!");
                }
                return(popms);
                break;
            }

        }
 ##       popms=c(popms,popm);
    }
}else{
    ampsm=ampsm_ob[which(ngroups_sm_ob==param$ngroup_show)];
    popm=plot_power(pcoa,rec,ampsm,param=param);
    print("Identified subpopulations:");
    print(ngroup);
    dev.set(3);
    points(ampsm,ngroup,pch=1,col="purple");
    dev.set(4);
    popms=append(list(popm),popms);
    names(popms)[1]=paste0("ngroup",ngroup);
    print(names(popms));
    return(popms);
}

    
}

##' @title popmap.out". 
#' @export
popmap.out <- function(ofile,amp_euc,amp_bin,omit_popid){
    if(length(omit_popid)>0){
        ofi=sprintf("%s_euc%d_bin%d_omit%s_",ofile,amp_euc,amp_bin,paste0(omit_popid,collapse="-"));
    }else{
        ofi=sprintf("%s_euc%d_bin%d_",ofile,amp_euc,amp_bin);
    }

        for(k in 2:5){
            pngout(k,plotfile=paste0(ofi,k),outeps=T,density=700,geometry=700);
        }

}


