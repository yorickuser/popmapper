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
    flag_calc_dist_all=1,
    flag_read_sample_location=1,
    flag_show_legend_location=1,
    flag_runid=1,
    flag_binary=2, ## 0 euclid, 1 my binary, 2 hybrid, 3 binary by R function
    nsamp=200,
    win_width=c(6,6,6,6),
    win_height=c(5,5,5,5),
    flag_pforeach=0,
    find_peak_scaling="maxmin",
    permutation_test=TRUE,
    edge_zmaxmin=0.95
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
find_peak <- function(pcoa,ampsm,ndiv=100,scaling="none"){
    px=pcoa$points[,1]
    py=pcoa$points[,2]
    ##    ampsm=0.08;


    
    sx=ampsm*(max(px)-min(px));
    ##sy=ampsm*(max(px)-min(px));
    xyratio=1.0;
    
    if(scaling=="maxmin"){
        xyratio=(max(py)-min(py))/(max(px)-min(px));
    }
    if(scaling=="mindis"){
        pxdis=px*0.0;
        for(i in 1:length(px))pxdis[i]=min(abs(px[-i]-px[i]));
        pxdis=mean(pxdis);
        pydis=py*0.0;
        for(i in 1:length(py))pydis[i]=min(abs(py[-i]-py[i]));
        pydis=mean(pydis);

        xyratio=pydis/pxdis;
    }

    sy=sx*xyratio;
    
    ##print(paste("scaling=",scaling," xyratio=",xyratio));
    
    
    
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
    cat("\nFraction of NULLs   mean:",p_NULL_mean,"   max:",max(p_NULL),"   min:",min(p_NULL), "\n");
    
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
        dist0=foreach::foreach(i=1:nn,.combine='rbind',.packages="foreach")%dopar%{    
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




calc_dist1 <- function(tab1,amp_euc=1,amp_bin=0,flag_binary=2,flag_pforeach=1,popids1=NULL,ids1=NULL){
    nn=nrow(tab1);
    nc=ncol(tab1);
    
    tab11=as.matrix(tab1);
    
    p_NULL=rowSums(tab11==0)/nc;
    p_NULL_mean=mean(p_NULL);
    cat("\n\nraction of NULLs   mean:",p_NULL_mean,"   max:",max(p_NULL),"   min:",min(p_NULL), "\n");
    
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
        dist0=foreach::foreach(i=1:nn,.combine='rbind',.packages="foreach")%dopar%{    
            calc_dist_unit(i,nn,tab11,mask,amp_euc,amp_bin,flag_binary=flag_binary);
        };
    }

    
    n_eu=dist0[,ncol(dist0)];
    dist0=dist0[,-ncol(dist0)];
    dist0=0.5*(dist0+t(dist0));

    nnu=nrow(dist0)%/%2;
    
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


    if(length(ids1)>0){
        ids2=ids1[seq(2,length(ids1),2)];
    }else{
        ids2=seq(nrow(dist0u));
    }
    if(length(popids1)>0){
        popids2=popids1[seq(2,length(popids1),2)];
    }else{
        popids2=rep(1,nrow(dist0u));
        }


    return(list(neu=n_eu,dist=dist0,dist2=dist0u,ids2=ids2,popids2=popids2));
}




  
read_data <- function(file_tsv=file_tsv,file_sample_location=file_sample_location,omit_popid=omit_popid,amp_euc=1,amp_bin=0,param=param0){
    if(flag_envs).ee.append("read_data",environment());

    if(length(file_sample_location)==0)param$flag_read_sample_location=0;
    
    print("reading file...");    
    df=read.csv(file_tsv,header=F,sep="\t",comment.char="#");
        
    tab0=df[-1,]
    tab0=as.matrix(tab0[,-c(1,2)]);
    tab0=apply(tab0,2,as.numeric);
    
    if(param$flag_random==1)tab0=tab0[,order(rnorm(ncol(tab0)))];
        
    sample_names=df[,1];
    sample_names=as.character(sample_names[-1]);
       
    if(param$flag_read_sample_location==1){
        sampleinfo=read.table(file_sample_location,header=TRUE,sep="\t");
        popids.si=sampleinfo$PopID;
        runids.si=sampleinfo$RunID;
        runlabels.si=as.character(sampleinfo$RunLabel);
        popnames.si=as.character(sampleinfo$PopName);
        
        ids.si=as.integer(gsub("S","",sampleinfo$SampleName))+as.integer(sampleinfo$RunID)*1000;

        runid=unique(runids.si);

        runlabel=rep("a",length(runid));

        for(i in 1:length(runid)){
            runlabel[i]=runlabels.si[(which(runids.si==runid[i]))[1]];
        }
        
        runlabel_buf=paste0(runlabel,"S");
        ##runid_with_s=paste0("S",seq(length(runlabel)));
        runid_with_s=paste0("S",runid);
        sample_namesb=sample_names;
        for(k in 1:length(runlabel_buf))sample_namesb=(gsub(runlabel_buf[k],runid_with_s[k],sample_namesb));
        ids=as.integer(gsub("S","",sample_namesb));

        
        ##if(flag_single_run==1)ids=ids+1000;
        if(max(ids)<1000)ids=ids+1000;
        
        
        ##SampleName2=paste0(runlabel[sampleinfo$RunID],sampleinfo$SampleName);
        SampleName2=paste0(sampleinfo$RunLabel,sampleinfo$SampleName);
        sampleinfo=data.frame(SampleName2=SampleName2,sampleinfo);
        
    }else{


        runlabels=substr(sample_names,1,1);

        sample_names.si=unique(sample_names);
        runlabels.si=substr(sample_names.si,1,1);
        
        runlabel=unique(runlabels.si);
        runid=seq(length(runlabel));

        
        runids.si=rep(0,length(sample_names.si));
        runids=rep(0,length(sample_names));
        for(i in 1:length(runid)){
            lis=which(runlabels.si==runlabel[i]);
            runids.si[lis]=runid[i];

            lis=which(runlabels==runlabel[i]);
            runids[lis]=runid[i];
        }


        ids=1000*runids+as.integer(substr(sample_names,3,100));
        
        ids.si=1000*runids.si+as.integer(substr(sample_names.si,3,100));
        
        ##        ids.si=unique(ids);
        ##ids.si=as.integer(gsub("S","",sampleinfo$SampleName))+as.integer(sampleinfo$RunID)*1000; 
        popids.si=rep(1,length(ids.si));
        ##runids=rep(1,length(ids.si));
        popnames.si=rep("pop1",length(ids.si));

        nsample=length(sample_names.si);
        sampleinfo=data.frame(SampleName2=sample_names.si,RunID=runids.si,RunLabel=runlabels.si,SampleName=substr(sample_names.si,2,100),PopName=rep("pop1",nsample),PopID=rep(1,nsample));

    }
    
    


    popids1=rep(0,nrow(tab0));
    runids1=popids1;
    mask0=rep(FALSE,nrow(tab0));
    
    for(i in 1:length(mask0)){
        lis=which(ids.si==ids[i]);
        if(length(lis)>0){
            mask0[i]=TRUE;
            popids1[i]=popids.si[lis];
            runids1[i]=runids.si[lis];
            }
    }

    if(sum(mask0)==length(mask0)){
        cat("All samples found their locations \n ");

    }
    if(sum(mask0)!=length(mask0)){
        cat("Samples lacking their locations are removed!!! :\n ");
        rlis=which(!mask0);
        cat(paste(rlis,sample_names[rlis],"\n"));        
    }
    
    tab1=tab0[mask0,];
    ids1=ids[mask0];
    sample_names1=sample_names[mask0];
    popids1=popids1[mask0];
    runids1=runids1[mask0];

    
    if(length(omit_popid)>0){
        omit_mask_popid=rep(FALSE,length(popids1));
        for(j in 1:length(omit_popid))omit_mask_popid=omit_mask_popid+(popids1==omit_popid[j]);
        
        omit_mask=(omit_mask_popid>0);
        
        ids_omit=ids1[omit_mask];
        popids1_omit=popids1[omit_mask];
        runids1_omit=runids1[omit_mask];
        sample_omit=sample_names1[omit_mask];
        
        cat("omit_popid: ",omit_popid,"\n");
        
        cat("omitted files:\n");
        print(paste(unique(sample_omit),collapse=", "));
        
        tab1=tab0[!omit_mask,];
        ids1=ids1[!omit_mask];
        popids1=popids1[!omit_mask];
        runids1=runids1[!omit_mask];
        sample_names1=sample_names1[!omit_mask];

    }

    omit=rep(0,nrow(sampleinfo));
    if(length(omit_popid)>0){
        for(j in 1:length(omit_popid))omit=omit+(popids.si==omit_popid[j]);
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

    return(list(tab00=tab00,tab1=tab1,tab1p=tab1p,ids1=ids1,popids1=popids1,runids1=runids1,popids.si=popids.si,popnames.si=popnames.si,ids.si=ids.si,runlabel=runlabel,runids.si=runids.si,runid=runid,sample_names1=sample_names1,sample.info=sampleinfo));
}


##' @title popmap.readdata". 
#' @export
popmap.readdata <- function(file_tsv=NULL,file_sample_location=NULL,  sample_group_name="Samples", omit_popid=NULL,amp_euc=1,amp_bin=0,param=param0){
    if(flag_envs).ee.append("popmap.readdata",environment());

    if(length(file_sample_location)==0)param$flag_read_sample_location=0;

  
    re_read_data=read_data(file_tsv=file_tsv,file_sample_location=file_sample_location, omit_popid=omit_popid,amp_euc=amp_euc,amp_bin=amp_bin,param=param);

        tab00=re_read_data$tab00;
        tab1=re_read_data$tab1;
        tab1p=re_read_data$tab1p;
        ids1=re_read_data$ids1;
    popids1=re_read_data$popids1;
    
        ids.si=re_read_data$ids.si;
        popids.si=re_read_data$popids.si;
    popnames.si=re_read_data$popnames.si;
    runlabel=re_read_data$runlabel;
    runid=re_read_data$runid;


    
    sample_names1=re_read_data$sample_names1;
    runids1=re_read_data$runids1;
        ##rm(re_read_data);
        
        
        print("calculating distance...");

    if(param$permutation_test==TRUE){
            res_calc_dist=calc_dist1(tab1,amp_euc=amp_euc,amp_bin=amp_bin,flag_binary=param$flag_binary,flag_pforeach=param$flag_pforeach,ids1=ids1,popids1=popids1);
            res_calc_distp=calc_dist1(tab1p,amp_euc=amp_euc,amp_bin=amp_bin,flag_binary=param$flag_binary,flag_pforeach=param$flag_pforeach,ids1=ids1,popids1=popids1);

            if(param$flag_calc_dist_all==1){
                res_calc_dist_all=calc_dist1(tab00,amp_euc=amp_euc,amp_bin=amp_bin,flag_binary=param$flag_binary,flag_pforeach=param$flag_pforeach,ids1=ids1,popids1=popids1);   
            }
            
    }else{
        res_calc_dist=calc_dist1(tab00,amp_euc=amp_euc,amp_bin=amp_bin,flag_binary=param$flag_binary,flag_pforeach=param$flag_pforeach,ids1=ids1,popids1=popids1);
        res_calc_distp=res_calc_dist;
        if(param$flag_calc_dist_all==1)res_calc_dist_all=res_calc_dist;
    }

        
    
        dist0=res_calc_dist$dist;
        dist2=res_calc_dist$dist2;
    n_eu=res_calc_dist$neu;
    ids2=res_calc_dist$ids2;
    popids2=res_calc_dist$popids2;
    
        
        dist0p=res_calc_distp$dist;
        dist2p=res_calc_distp$dist2;
        n_eup=res_calc_distp$neu;

    dist2.all=NULL;
    dist0.all=NULL;
    if(param$flag_calc_dist_all){
        dist2.all=res_calc_dist_all$dist2;
        dist0.all=res_calc_dist_all$dist;
}
    ##nrun=as.integer(max(ids2)/1000);
    nrun=max(runid);
    masks=matrix(ids2 > -1, nrow=nrun, ncol=length(ids2), byrow=T);
    for(i in 1:nrun)masks[i,]=(i*1000<=ids2)+(ids2<(i+1)*1000)==2;
    mask=colSums(masks)>0;

    
    
        popid_u=unique(popids.si);
        popname_u=rep("",length(popid_u));
        for(i in 1:length(popid_u))popname_u[i]=popnames.si[(which(popids.si==popid_u[i]))[1]];
        
        popnames2=popname_u[popids2];

        sample_names2=rep("a",length(ids2));
        for(i in 1:length(ids2)){
            lis=which(ids1==ids2[i]);
            sample_names2[i]=sample_names1[lis[1]];
        }

    lis2=ids2*0;
        for(i in 1:length(ids2)){
        lis2[i]=(which(ids.si==ids2[i]))[1];
    }
    sample.info2=re_read_data$sample.info[lis2,];

    
        return(list(dist2=dist2,dist2p=dist2p,ids2=ids2,popids2=popids2,popnames2=popnames2,sample_names2=sample_names2,sample_names1=sample_names1,popid.legend=popid_u,popname.legend=popname_u,tab00=tab00,tab1=tab1,tab1p=tab1p,ids1=ids1,popids1=popids1,runids1=runids1,ids.si=ids.si,popids.si=popids.si,popnames.si=popnames.si,neu=n_eu,neup=n_eup,masks=masks,mask=mask,nrun=nrun, runlabel=runlabel,runid=runid,sample.info=re_read_data$sample.info,sample.info2=sample.info2,amp_euc=amp_euc,amp_bin=amp_bin,omit_popid=omit_popid,sample_group_name=sample_group_name,dist2.all=dist2.all,dist0.all=dist0.all));   
}

##get_pcoa2d <- function(dist2,nsamp,ampst=0.8,amped=1.2){
get_pcoa2d <- function(dist2,nsamp,ampst=0.1,amped=1.2){
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

    
    return(list(points=res$points,dis2=dis2,p1=p1,p2=p2,edges=edges,eigenvalues=eig));
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
    ngroups2=foreach::foreach(i=1:nn,.combine='rbind',.packages="foreach")%dopar%{
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

calc_groups_sm <- function(dist2,edges,edge2ampsm,ndiv=100,flag_pforeach=param0$flag_pforeach,scaling=param0$find_peak_scaling){

    if(flag_envs).ee.append("calc_groups_sm",environment());

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
            (find_peak(res,edges[i]*edge2ampsm,ndiv=ndiv,scaling=scaling))$ngroup; 
        });
    }else{
        print("using foreach..");
        ngroups_sm=foreach::foreach(i=1:nn,.combine='rbind',.packages="foreach")%dopar%{
            if(i%%50==0)cat(100*i/nn,"% ");
            (find_peak(res,edges[i]*edge2ampsm,ndiv=ndiv,scaling=scaling))$ngroup;
            
        };
    }

    
    print("done");
    
    ngroups_sm_ob=sort(unique(ngroups_sm));
    ampsm_ob=ngroups_sm_ob*0.0;
    for(i in 1:length(ngroups_sm_ob)){

        edgebufs= edges[ngroups_sm==ngroups_sm_ob[i]];
        if(length(edgebufs)>1){
            ##ampsm_ob[i]=edge2ampsm*edgebufs[as.integer(length(edgebufs)%/%2)];
            ampsm_ob[i]=edge2ampsm*min(edgebufs);
        }else{
            ampsm_ob[i]=edgebufs;
        }
##        ampsm_ob[i]=edge2ampsm*min(edges[ngroups_sm==ngroups_sm_ob[i]]);
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

map_power <- function(pcoa,rec,ampsm_opt,ndiv=100,scaling=param0$find_peak_scaling){
    if(flag_envs).ee.append("map_power",environment());

    popids2=rec$popids2;
   
    px = pcoa$points[,1]
    py = pcoa$points[,2]
   
    ampsm = ampsm_opt;
    
    res_peak=find_peak(pcoa,ampsm_opt,ndiv=ndiv,scaling=scaling);
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

    zmaxmin=matrix(rep(0.0,3*ngroup),nrow=3);
    for(k in 1:ngroup){
        emask=(ccer1[2:(ndiv-1),2:(ndiv-1),k]==1)*(ccer1[1:(ndiv-2),2:(ndiv-1),k]+ccer1[3:(ndiv),2:(ndiv-1),k]+ccer1[2:(ndiv-1),1:(ndiv-2),k]+ccer1[2:(ndiv-1),3:(ndiv),k]!=4);

        zmaxmin[1,k]=max(zz[ccer1[,,k]==1]);
        zmaxmin[2,k]=max((zz[2:(ndiv-1),2:(ndiv-1)])[emask==1]);
        zmaxmin[3,k]=zmaxmin[2,k]/zmaxmin[1,k];
        
    }
    
##    print("AAA");

  ##  print(lispeak);
    return(list(cid_sm=cid_sm,peakx=peakx,peaky=peaky,xx0=xx0,yy0=yy0,xx=xx,yy=yy,zz=zz,sx=sx,sy=sy,ccer0=ccer0,ccer1=ccer1,ndiv=ndiv,zmaxmin=zmaxmin));
}


pfunc <- function(px,py,xx0,yy0,peakx,peaky,ccer0,ccer1,cid_sm,labels,lev=0.95,pals= NULL,flag_map=1,mask=NULL,col_mask="white",peak=TRUE,power=TRUE,pals_power=NULL){

    if(flag_envs).ee.append("pfunc",environment());

  
    if(flag_map==0){
        text(x=px,y=py,col=pals[cid_sm],labels=labels,cex=0.7,font=1)
    }else{

        
        if(power==TRUE){
            for(ii in 1:length(pals_power)){
                ##        points(xx[which(ccer1[,,ii]>0)],yy[which(ccer1[,,ii]>0)],col=pals[ii],pch=20,cex=0.7);
                contour(xx0,yy0,ccer1[,,ii],add=TRUE,drawlabels=FALSE,method="simple",col=pals_power[ii],levels=c(lev),labels=NULL,lwd=2,lty=1)
                contour(xx0,yy0,ccer0[,,ii],add=TRUE,drawlabels=FALSE,method="simple",col=pals_power[ii],levels=c(lev),labels=NULL,lwd=1,lty=2)   
                ##        text(xx[which(ccer1[,,ii]>0)],yy[which(ccer1[,,ii]>0)],label=paste(ii),col=pals[ii],pch=20,cex=0.5);
            }
        }



        text(x=px,y=py,col=pals[cid_sm],labels=labels,cex=0.7,font=1)
        if(length(mask)==0){
            
            if(peak){
                points(peakx,peaky,col="black",pch=17,lwd=3,cex=1.8);
                points(peakx,peaky,col=pals,pch=17,lwd=3,cex=1.5);
                points(peakx,peaky,col="white",pch=17,lwd=3,cex=0.8);
            }
        }else{
            lis0=(cid_sm!=mask);
            lis1=(cid_sm==mask);
           ## text(x=px,y=py,col="white",labels=labels,cex=0.7,font=1)
           ## text(x=px[lis1],y=py[lis1],col=pals[cid_sm[lis1]],labels=labels[lis1],cex=0.7,font=1);
            if(peak){
            points(peakx[mask],peaky[mask],col="black",pch=1,lwd=3,cex=1.8);
            points(peakx[mask],peaky[mask],col="white",pch=1,lwd=3,cex=1.2);
            points(peakx[mask],peaky[mask],col=pals[mask],pch=1,lwd=3,cex=0.6);
            }

        }
        
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
    if(flag_envs).ee.append("cmdist",environment());
    
    ngroup=max(cid_sm);
    gdis=matrix(rep(0.0,ngroup*ngroup),nrow=ngroup,ncol=ngroup);
    pvs=gdis;
    dij=gdis;
    group_size=rep(0,ngroup);
    for(i in 1:ngroup)group_size[i]=sum(cid_sm==i);
    minsizes=sort(group_size)[1:2];
    print("##group sizes:");
    print(seq(ngroup));
    print(group_size);
    
    if(sum(minsizes)<=3){
        test=FALSE;
        print("!!There are two small groups!!");
    }
    
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
             cmdist_unit(k, ilis,jlis,dist2p,cid_sm,test=test,flag_welch_t_test=flag_welch_t_test,nperm=nperm,flag_pforeach=flag_pforeach);
         });
         
     }else{
         print("using foreach..");
         res=foreach::foreach(k=1:length(ilis),.combine='rbind',.packages="foreach")%dopar%{
             cmdist_unit(k, ilis,jlis,dist2p,cid_sm,test=test,flag_welch_t_test=flag_welch_t_test,nperm=nperm,flag_pforeach=flag_pforeach);
         };

     }
         
##    print(c(length(ilis),res));

    for(k in 1:length(ilis)){
        if(length(ilis)==1){
            dm0=res[1];
            pv0=res[2];
            dij0=res[3];
            
        }else{
            dm0=res[k,1];
            pv0=res[k,2];
            dij0=res[k,3];
            
        }
        
        gdis[ilis[k],jlis[k]]=dm0;
        gdis[jlis[k],ilis[k]]=dm0;
        pvs[ilis[k],jlis[k]]=pv0;
        pvs[jlis[k],ilis[k]]=pv0;
        dij[ilis[k],jlis[k]]=dij0;
        dij[jlis[k],ilis[k]]=dij0;
  
    }
    ##}
    

    ##    print("a");
    ##    print(pvs);
    return(list(gdis=gdis,pvs=pvs,dij=dij));
}


#' @export
cmdist_unit <- function(k, ilis,jlis,dist2p,cid_sm,test=FALSE,flag_welch_t_test=0,nperm=1000,flag_pforeach){
##       if(flag_envs).ee.append("cmdist_unit",environment());
    i=ilis[k];
    j=jlis[k];
    lis=which((cid_sm==i)+(cid_sm==j)>0);
    dist2pp=dist2p[lis,lis];
    cid=cid_sm[lis];
    rec=cmdscale(dist2pp);
    
       
       lisi=which(cid==i);
       lisj=which(cid==j);
       nni=length(lisi);
       nnj=length(lisj);

        cat("##",nni,":",nnj,"*gid*",unique(cid));
        ppi=rec[lisi,,drop=F];
        ppj=rec[lisj,,drop=F];
        ppim=apply(ppi,2,mean);
        ppjm=apply(ppj,2,mean);
        
        covi=0.0;
        covj=0.0;
        if(nni>1){
            covi=cov(ppi);
     
        }
        if(nnj>1){
            covj=cov(ppj);
          
        
        }
        

        ma=(covi*(nni-1)+covj*(nnj-1))/(nni+nnj-2);

        if(nni+nnj>3){
            ev=solve(ma)%*%(ppim-ppjm);
            
        }else{
            ev=ppim-ppjm;
        }
        
       ev=ev/sqrt(ev[1]^2+ev[2]^2);
       ev0=ev;
       ppi0=ppi;
       ppj0=ppj;
       ppim0=ppim;
       ppjm0=ppjm;
       
        cat(sprintf("\n gid %d:%d direction adjust for between-group distances ev=(%2.3f,%2.3f)\n",i,j,ev[1],ev[2]));
        
        pp=rec[,1]*ev[1]+rec[,2]*ev[2];      
        pi=ppi[,1]*ev[1]+ppi[,2]*ev[2];
        pj=ppj[,1]*ev[1]+ppj[,2]*ev[2];
    
    
    
    sdi=sd(pi);
    sdj=sd(pj);
    dij0=abs(mean(pi)-mean(pj));
    pv=1.0;
    ##sdw=(length(pi)*sdi+length(pj)*sdj)/(length(lis));
    sdw0=sqrt((sqm(pi)+sqm(pj))/(nni+nnj-2));
    
    if(flag_welch_t_test==1)sdw1=sqrt(sdi^2/length(pi)+sdj^2/length(pj));
    
    ##gdis[i,j]=dij/(0.5*(sdi+sdj));
    ##gdis[i,j]=dij/sdw;
    ##gdis[j,i]=gdis[i,j];

       dm0=dij0/sdw0;

      use_pforeach=flag_pforeach;
    
    if(test==TRUE){
        if(min(nni,nnj)>1){
            if(flag_welch_t_test==1)sdw0=sdw1;
                inde0=dij0/sdw0;
                
                
                ##nperm=1000;
            indes=rep(0.0,nperm);
            indes0=rep(0.0,nperm);
            rns=Reshape((rnorm(length(cid)*nperm)),length(cid),nperm);
            
            neach=nperm%/%10;
            ##neach=1000;
            nk=nperm%/%neach;
            ke=c(seq(nk)*neach,(nk*neach+nperm%%neach));
            ks=c(1,seq(nk)*neach+1);
            nkk=nk+as.integer(nperm%%neach>0);

            if(use_pforeach==1){
                print("pforeach for each test");
                ##                for(kk in 1:nperm){
                ##indes=pforeach::pforeach(kk=1:nperm, .combine='c',.errorhandling="stop")({
                
                indes=pforeach::pforeach(kk=1:nkk, .combine='c',.errorhandling="stop")({
                    klis=(ks[kk]:ke[kk]);
                    rns1=rns[,klis];
                    npf=length(klis);
                    indes1=rep(0.0,npf);
                    
                    for(m in 1:npf){                        
                        lisr=order(rns1[,m]);                      
                        indes1[m]=popmapper::get_inde(i,j,cid[lisr],rec,nni,nnj,flag_welch_t_test=flag_welch_t_test);
                    }
                    indes1;
                
                });

                
            }else{
                print("foreach for each test");
                

                indes=foreach::foreach(kk=1:nkk,.combine='c',.packages="foreach")%dopar%{
                    klis=(ks[kk]:ke[kk]);
                    rns1=rns[,klis];
                    npf=length(klis);
                    indes1=rep(0.0,npf);
                    
                    for(m in 1:npf){                        
                        lisr=order(rns1[,m]);                      
                        indes1[m]=popmapper::get_inde(i,j,cid[lisr],rec,nni,nnj,flag_welch_t_test=flag_welch_t_test);
                    }
                    indes1;
                
                };
                
                
               

            }
            
            ##        indes[kk]=dij/sdw;
            ##    }
                p_value=(sum(indes>inde0)+sum(indes==inde0)*0.5)/nperm;

            cat(sprintf("pvalue:%f  distance:%f \n", p_value, dm0));
            ##cat("!!!!!");
            ##print(c(sd(indes),max(indes),min(indes),inde0,p_value));
            ##print(indes-indes0);
            
                
                
                if(!is.na(p_value)){
                     pv=p_value;
                     ##pvs[j,i]=p_value;
                }
        }
        }
    dij00=abs(mean(rec[lisi,1])-mean(rec[lisj,1]));
    ##dij00=abs(mean(dist2pp[lisi,lisj]));
            return(c(dm0,pv,dij00))
     }

#' @export
get_inde <- function(i,j,cid1,rec,nni,nnj,flag_welch_t_test=0){
    lisi=which(cid1==i);
    lisj=which(cid1==j);
    
    ppi=rec[lisi,];
    ppj=rec[lisj,];
    
    ppi=rec[lisi,];
    ppj=rec[lisj,];
    
    ma=(cov(ppi)*(nni-1)+cov(ppj)*(nnj-1))/(nni+nnj-2);
    
    ev=solve(ma)%*%(apply(ppi,2,mean)-apply(ppj,2,mean));
    ev=ev/sqrt(ev[1]^2+ev[2]^2);
    
    pi=ppi[,1]*ev[1]+ppi[,2]*ev[2];
    pj=ppj[,1]*ev[1]+ppj[,2]*ev[2];
                            
                            sdi=sd(pi);
                            sdj=sd(pj);
                            dij=abs(mean(pi)-mean(pj));
                            sdw=sqrt((sqm(pi)+sqm(pj))/(nni+nnj-2));
                           if(flag_welch_t_test==1){
                                sdw1=sqrt(sdi^2/nni+sdj^2/nnj);
                                sdw=sdw1;
                            }
                            inde1=dij/sdw;
          
                            return (inde1);
   }

##' @title plot.popmap". 
#' @export
plot.popmap <- function(popmap,param=param0,perm=0,resp=NULL,main=NULL,labels=FALSE){
    popmap.plot(popmap,param=param,perm=perm,resp=resp,main=main,labels=labels);
}

##' @title popmap.plot". 
#' @export
popmap.plot <- function(popmap,param=param0,perm=0,resp=NULL,main=NULL,labels=NULL,sample_group_name="",peak=TRUE,power=TRUE,main_add=FALSE,cid_sm=popmap$cid_sm,pals=popmap$pals){
    if(flag_envs).ee.append("popmap.plot",environment());
    if(main_add)main0=main;

    cid_sm_power=cid_sm;
    
    if(length(pals)!=length(cid_sm)){
        cid_sm_power=popmap$cid_sm;    
        pals_power=pals;
        if(length(pals)!=max(cid_sm))pals=rainbow(max(cid_sm),v=1.0);   
        ##print(pals);
        ##print(pals_power);
    }
    ##cid_sm=popmap$cid_sm;
    if(length(labels)==1){
        if(class(labels)=="character"){
            lcase=c("gid","gid.merge","gid.sub","popid");
            labels0=list(gid=paste(popmap$cid_sm),gid.merge=paste(popmap$cid_sm.merge),gid.merge=paste(popmap$cid_sm.sub),popid=popmap$popids2);
            labels=labels0[[which(lcase==labels)]];
        }
    }
    if(length(labels)==0){
        labels=popmap$popids2;
        
    }
    if(perm==1){
        pcoap=popmap$pcoap;
        powermap=popmap$powermap;

        if(length(resp)==0)resp=find_peak(pcoap,popmap$ampsm,ndiv=powermap$ndiv,scaling=param$find_peak_scaling);

        if((length(main)==0)||(main_add)){
            if(popmap$maxp>0){
                maxij=(which(popmap$gdist$pvs==max(popmap$gdist$pvs),arr.ind=T))[1,];
                main=sprintf("%s  max p-value: %f (%d:%d nperm: %d)",sample_group_name,max(popmap$gdist$pvs),maxij[1],maxij[2],param$nperm);
                if(main_add)main=paste(main0,main);
                
            }else{
                main=sprintf("%s  max p-value: %f (nperm: %d)",sample_group_name,popmap$maxp,param$nperm);
             
            }
            
        }

        
        plot_power_map(pcoap,resp,cid_sm,flag_fix_lim=param$flag_fix_lim,pals=pals,ccer0=powermap$cerr0,ccer1=powermap$ccer1,labels=labels,main=main,flag_map=0,param=param,peak=peak,power=power,cid_sm_power=cid_sm_power,pals_power=pals_power);

        
    }else{
        pcoa=popmap$pcoa;
        powermap=popmap$powermap;
        if((length(main)==0)||(main_add))main=sprintf("%s  ngroups: %d   scaling: %s  sd: %2.3f\n  max p-value: %f    nperm: %d  zmaxmin: %2.3f",sample_group_name,popmap$ngroup, param$find_peak_scaling,popmap$ampsm,popmap$maxp,param$nperm, max(popmap$zmaxmin[3,]));
        if(main_add)main=paste(main0,main);


        plot_power_map(pcoa,powermap,cid_sm,flag_fix_lim=0,pals=pals,ccer0=powermap$ccer0,ccer1=powermap$ccer1,labels=labels,flag_map=1,main=main,param=param,peak=peak,power=power,cid_sm_power=cid_sm_power,pals_power=pals_power);
        
        
    }
    invisible(pals);
}

plot_power_map <- function(pcoa,res,cid_sm,flag_fix_lim=0,pals=NULL,ccer0=NULL,ccer1=NULL,labels=NULL,flag_map=1,main=NULL,param=param0,mask=NULL,col_mask="white",peak=TRUE,power=TRUE,cid_sm_power=cid_sm,pals_power=NULL){
    if(flag_envs).ee.append("plot_power_map",environment());
    xx0=res$xx0;
    yy0=res$yy0;
    xx=res$xx;
    yy=res$yy;
    zz=res$zz;
    sx=res$sx;
    sy=res$sy;


##    zmm=max(zmaxmin[3,]);
    
    ngroup = max(cid_sm);
    peakx = res$peakx;
    peaky = res$peaky;
    
    px=pcoa$points[,1];
    py=pcoa$points[,2];

    ngroup=max(cid_sm);
    
    if(max(cid_sm)!=max(cid_sm_power)){
      
             for(i in 1:length(pals_power)){
                 lis=which(cid_sm_power==i);
                 mid=(length(lis)+1)%/%2;
                 pals1=unique(pals[cid_sm[lis]]);
                 mid=(length(pals1)+1)%/%2;
                 
                pals_power[i]=pals1[mid];
                ##pals_power[i]=pals[cid_sm[lis[1]]];
                ##pals_power[i]=pals[i];
                
            }
        }
        
   
        
    par(mar=c(5,5,5,6));
    par(xpd=T);
    dev.hold();
    
    if(length(pals)==0)pals = rainbow(ngroup,v=1.0);
    if(length(pals_power)==0)pals_power=pals;
    
    nlev=20;
    
    if(param$flag_fix_lim==0){
        plot(xx0,yy0,type="n",xlim=c(min(xx0),max(xx0)),ylim=c(min(yy0),max(yy0)),main=main,xaxs="i",yaxs="i",xlab="Axis 1", ylab="Axis 2",cex.main=1.1);
    }else{    
        plot(xx0,yy0,type="n",xlim=xlim0,ylim=ylim0,main=main,xaxs="i",yaxs="i",xlab="Axis 1", ylab="Axis 2",cex.main=1.1);
        
        polygon(c(xlim0[1],xlim0[2],xlim0[2],xlim0[1]),c(ylim0[1],ylim0[1],ylim0[2],ylim0[2]), col=(gray.colors(nlev,start=0.1,end=0.6))[1]);
        
    }
    

    .filled.contour(xx0,yy0,zz,levels=seq(0,max(zz),,nlev),col=gray.colors(nlev,start=0.1,end=0.6));


            
    pfunc(px,py,xx0,yy0,peakx,peaky,ccer0,ccer1,cid_sm,labels,pals=pals,flag_map=flag_map,mask=mask,col_mask=col_mask,peak=peak,power=power,pals_power=pals_power);
    
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

plot_power <- function(pcoa,rec,ampsm_opt,param=param0,sample_group_name=""){
    if(flag_envs).ee.append("plot_power",environment());
    dist2=rec$dist2;
    dist2p=rec$dist2p;
    dist2.all=rec$dist2.all;
    popids2=rec$popids2;

    ndiv=param$ndiv;
    
    res=map_power(pcoa,rec,ampsm_opt,ndiv=param$ndiv,scaling=param$find_peak_scaling);
    cid_sm = res$cid_sm;
    ngroup = max(cid_sm);

    zmaxmin=res$zmaxmin;
    
    dev.set(4);

    main=sprintf("%s sd: %2.3f     groups: %d   zmaxmin: %2.3f",sample_group_name,ampsm_opt,ngroup,max(zmaxmin[3,]));

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

##        edge_zmaxmin=0.9;
        if(max(zmaxmin[3,])>param$edge_zmaxmin){
            print("Too shallow valley!!!");
            print(zmaxmin);
            flag_do_perm=0;
        }
        gdist=NULL;
        gdist0=NULL;
        gdist.all=NULL;
        if(flag_do_perm==1){
            print("calculating genetic distances among subpopulations..");
            gdist=cmdist(dist2p,cid_sm,test=param$permutation_test,flag_welch_t_test=param$flag_welch_t_test,nperm=param$nperm,flag_pforeach=param$flag_pforeach);
            
            gdist0=cmdist(dist2,cid_sm,test=FALSE,flag_welch_t_test=param$flag_welch_t_test,nperm=param$nperm,flag_pforeach=param$flag_pforeach);


            if(length(dist2.all)>0)gdist.all=cmdist(dist2.all,cid_sm,test=FALSE,flag_welch_t_test=param$flag_welch_t_test,nperm=param$nperm,flag_pforeach=param$flag_pforeach);
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
        
        
        resp=find_peak(pcoap,ampsm_opt,ndiv=ndiv,scaling=param$find_peak_scaling);
        

        maxij=1.0;
        if(flag_do_perm==1){
            maxij=(which(gdist$pvs==max(gdist$pvs),arr.ind=T))[1,];
        }    
        if(flag_do_perm==1){
            main0=sprintf("max p-value: %f (%d:%d nperm: %d)",max(gdist$pvs),maxij[1],maxij[2],param$nperm);
        }else{
            main0="too small populations!!";
        }

        
        ##plot_power_map(pcoap,resp,cid_sm,flag_fix_lim=param$flag_fix_lim,pals=pals,ccer0=res$cerr0,ccer1=res$ccer1,labels=paste(cid_sm),main=main0,flag_map=0,param=param);
        plot_power_map(pcoap,resp,cid_sm,flag_fix_lim=param$flag_fix_lim,pals=pals,ccer0=res$cerr0,ccer1=res$ccer1,labels=popids2,main=main0,flag_map=0,param=param);

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
    
    return(list(ngroup=ngroup,maxp=maxp,cid_sm=cid_sm,popids2=popids2,gdist=gdist,gdist0=gdist0,gdist.all=gdist.all,powermap=res,pcoa=pcoa,pcoap=pcoap,ampsm=ampsm_opt,pals=pals,zmaxmin=zmaxmin));
}


pngout<-function(dev_id=as.numeric(dev.list()[length(dev.list())]),plotfile="testout",density=150,geometry=600,outeps=FALSE,prefix="./",show=TRUE,outpng=TRUE,alpha=FALSE){
    
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
        if(alpha==FALSE){
            system(sprintf("convert -density %sx%s -geometry %s  -background white -alpha remove %s.eps %s.png",dens,dens,geom,tmpf,tmpf));
        }else{
            system(sprintf("convert -density %sx%s -geometry %s   %s.eps %s.png",dens,dens,geom,tmpf,tmpf));
        }
        system(sprintf("mv %s.png %s%s.png",tmpf,prefix,plotfile));
        cat(sprintf("png output: %s%s.png \n",prefix,plotfile));
        if(show)system(sprintf("display %s%s.png&",prefix,plotfile));
    
    }
    system(sprintf("rm %s.eps",tmpf));
    
}


##' @title popmap.plotdata". 
#' @export
popmap.plotdata <- function(pcoa,data,param=param0,sample_group_name=NULL,flag_lim_out=0,plot_mask=NULL,col1="red",col2="gray",tagid=NULL){
   if(flag_envs).ee.append("popmap.plotdata",environment());

   if(length(sample_group_name)==0)sample_group_name=data$sample_group_name;

   amp_euc=data$amp_euc;
   amp_bin=data$amp_bin;
   omit_popid=data$omit_popid;
   
   runlabel=data$runlabel;
   masks=data$masks;

   if(length(plot_mask)==0){
       mask=data$mask;
   }else{
       mask=plot_mask;
   }


    nrun=data$nrun;
    popnames.si=data$popnames.si;
    popids2=data$popids2;
    popids.si=data$popids.si;
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

    if(length(plot_mask)==0){
        text(x=pcoa$points[mask,1],y=pcoa$points[mask,2],col=cols[mask],labels=popids2,cex=0.7,font=1);


    }else{
        text(x=pcoa$points[mask,1],y=pcoa$points[mask,2],col=col1,labels=popids2[mask],cex=0.7,font=1);
        text(x=pcoa$points[!mask,1],y=pcoa$points[!mask,2],col=col2,labels=popids2[!mask],cex=0.7,font=1);
        }
    

        run_name="run";
        if(param$flag_cols_cerwave==1)run_name="group";
        if(length(plot_mask)==0){
            legend(param$legend_side, legend = paste(run_name,data$runlabel), col = pal,pch=rep(15,length(pal)),cex=0.8);
        }else{
            legend(param$legend_side, legend = tagid, col = col1,pch=15,cex=0.8);

            }
    
##    legend(par()$usr[2], par()$usr[3], legend = paste("run",seq(length(pal))), col = pal,pch=rep(15,length(pal)))

if(param$flag_show_legend_location==1){
##if(param$flag_read_sample_location==1){
    ##popid_u=unique(popids);
    ##popname_u=rep("",length(popid_u));
    ##for(i in 1:length(popid_u))popname_u[i]=popnames[(which(popids==popid_u[i]))[1]];
    
    
    legend(par()$usr[2], par()$usr[4], legend = paste0(data$popid.legend,": ",data$popname.legend), col = "black",cex=0.7,bty="n");

    text(par()$usr[2], par()$usr[3], labels = paste("Omitted \npopids:",paste(omit_popid,collapse=", ")), col = "red",cex=0.8,pos=4)

    ##mtext(paste0("Omitted popids:",omit_popid), side=4,col = "red",cex=0.8)
  ##  }
}
}

if(param$flag_runid==0){
    plot.new();    
    plot(pcoa$points,type="n",main="run_id: number,   population_id: color",xlab="Axis 1");

    pops=sort(unique(popids.si));
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
popmap.find <- function(rec,param=param0,sample_group_name=NULL,plot_data=1){
    if(flag_envs).ee.append("popmap.find",environment());

    if(length(rec$sample_names2)==0){
        rec$sample_names2=sprintf("S%03d",seq(dim(data$dist2)[1]));
    }

    if(length(rec$popids2)==0){
        rec$popids2=seq(length(rec$sample_names2));
        }
   if(length(rec$sample.info)==0){
                    rec$sample.info=data.frame(SampleName2=rec$sample_names2,PopID=rec$popids2);
                }  
    
    if(length(sample_group_name)==0)sample_group_name=rec$sample_group_name;
    if(length(dev.list())==0)xinit(param$win_width,param$win_height);
    
    dist2=rec$dist2;
    dist2p=rec$dist2p;    
    pcoa=get_pcoa2d(dist2,param$nsamp);
    edges=pcoa$edges;
    
    print("find groups for different edge values..");
    ngroups2=calc_groups(dist2,edges,param=param);
    
    print("find groups (peaks) for different gaussian smoothing..");
    re_calc_groups_sm=calc_groups_sm(dist2,edges,param$edge2ampsm,ndiv=param$ndiv,flag_pforeach=param$flag_pforeach,scaling=param$find_peak_scaling);
    ngroups_sm=re_calc_groups_sm$ngroups_sm;
    ampsm_ob=re_calc_groups_sm$ampsm_ob;
    ngroups_sm_ob=re_calc_groups_sm$ngroups_sm_ob;

    
############################################################################
    if(plot_data==1){
        dev.set(2);
        
        lims=popmap.plotdata(pcoa,rec,param=param,sample_group_name=sample_group_name,flag_lim_out=1);
        xlim0=lims$xlim0;
        ylim0=lims$ylim0;
    }
#####################################################################

dev.set(3);
    edge2ampsm=param$edge2ampsm;
    edgebuf=c(edges,edges*edge2ampsm);
plot(edges*edge2ampsm,ngroups_sm,type="l",log="xy",col="blue",ylim=c(1,ncol(dist2)),xlim=c(min(edgebuf),max(edgebuf)),xlab="sd for Gaussian smoothing",ylab="Number of groups")
lines(edges,ngroups2,col="black");

    if(param$permutation_test==TRUE){
        legend("topright", legend = c("Gaussian smooth","Neighbor connect","Significance"), col = c("blue","black","red"),pch=c(15,15,1),cex=0.8);
    }else{
        legend("topright", legend = c("Gaussian smooth","Neighbor connect","Obtained Popmap"), col = c("blue","black","red"),pch=c(15,15,1),cex=0.8);
    }
    
    popms=list(pcoa=pcoa);
    if(param$ngroup_show<0){
        pvs=rep(1.0,20);
        
        
        count_sig=0;
        for(k in 2:20){
            cat("\n\n::::Testing ",k, "groups-mapping  (ampsm:",ampsm_ob[k],") \n");
            maxp=1.0;
            popm=plot_power(pcoa,rec,ampsm_ob[k],param=param,sample_group_name=sample_group_name);
            maxpb=popm$maxp;
            
            if(!is.na(maxpb))maxp=maxpb;
            pvs[k]=maxp;
            
            if((maxp<param$edge_p_value)||(param$permutation_test==FALSE)){
                count_sig=count_sig+1;
                ##            si=data.frame(rec$sample.info,GroupID=rep(0,nrow(rec$sample.info)),pcoa.x=pcoa$points[,1],pcoa.y=pcoa$points[,2]);
             
                    buf=rep(0,nrow(rec$sample.info));
                    si=data.frame(rec$sample.info,GroupID=buf,pcoa.x=buf,pcoa.y=buf);

                sname=rec$sample_names2;
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
                class(popm) = 'popmap';
                popms=append(list(popm),popms);
                names(popms)[1]=paste0("ngroup",ngroup);
            
                dev.set(3);
            
                if(param$permutation_test){
                    points(ampsm_ob[k],ngroups_sm_ob[k],pch=1,col="red",lwd=1);
                }else{
                    points(ampsm_ob[k],ngroups_sm_ob[k],pch=1,col="red",lwd=1);
                }
                
                dev.set(4);

                if((ngroups_sm_ob[k]>=param$group_max)||(k==length(ngroups_sm_ob))){
                    print(names(popms));
                    return(popms);
                    break;
                    }
            }else{
                if(param$permutation_test==TRUE)print(sprintf("p-value %2.3f higher than %2.3f  ->>  not significant !!!",maxp,param$edge_p_value));
                
                
                if((ngroups_sm_ob[k]>param$group_max)||(k==length(ngroups_sm_ob))){
                    
                        if(count_sig>0){
                        ##opid=max(which(pvs<param$edge_p_value));
                        ##popm=plot_power(pcoa,rec,ampsm_ob[opid],param=param);
                        
                        dev.set(4);
                        popmap.plot(popms[[1]],param=param,perm=0,sample_group_name=sample_group_name)
                        dev.set(5);
                        popmap.plot(popms[[1]],param=param,perm=1)
                        
                        ##maxpb=popm$maxpb;
                        cat("\n\n::::Maximum identified subpopulations: ",popms[[1]]$ngroup,"\n");
                        if(param$permutation_test==TRUE)cat("Max p-value: ",popms[[1]]$maxp,"\n");
                        
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
    ngroup=popm$ngroup;
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



#####################################################################
#####################################################################

subdiv_edge_size <<- 4;
#' @export
popmap.out.gdis <- function(map,ofile,amp_euc=1,amp_bin=1,omit_popid=NULL){
    if(length(omit_popid)>0){
        ofi=sprintf("%s_euc%d_bin%d_omit%s_",ofile,amp_euc,amp_bin,paste0(omit_popid,collapse="-"));
    }else{
        ofi=sprintf("%s_euc%d_bin%d_",ofile,amp_euc,amp_bin);
    }

    if(length(map$gdis.all)>0){
        gdis=map$gdis.all;
    }else{
        print("gdis.all is approximately calculated as 0.5*(map$gdist$gdis+map$gdis0$gdis");
        gdis=0.5*(map$gdist$gdis+map$gdist0$gdis);
    }
    
    ofi_gdis=paste0(ofi,"gdis.tsv");
    ofi_pvs=paste0(ofi,"pvs.tsv");
    ofi_sampleinfo=paste0(ofi,"sample.info.tsv");
    
    write.table(map$sample.info,ofi_sampleinfo,sep="\t",col.names=T,row.names=F,quote=F);
    cat("popmap:\n",ofi_sampleinfo,"\n");
    
    write.table(gdis,ofi_gdis,sep="\t",col.names=F,row.names=F,quote=F);
    cat("genetic distances:\n",ofi_gdis,"\n");
    write.table(map$gdist$pvs,ofi_pvs,sep="\t",col.names=F,row.names=F,quote=F);
    cat("pvalues:\n",ofi_pvs,"\n");

    if(length(map$gdist.all.sub)>0){
        gdis=map$gdist.all.sub$gdis;
        ofi_gdis_sub=paste0(ofi,"gdis_sub.tsv");
        
        ##write.csv(round(gdis,digits=3),ofi_gdis_sub);
        ##write("\t",file=ofi_gdis_sub);
        colnames(gdis)[1]=paste0("\t",colnames(gdis)[1]);
        write.table(round(gdis,digits=3),file=ofi_gdis_sub,sep="\t",row.names=T,col.names=T,quote=F);
        cat("genetic distances among all subdivided groups:\n",ofi_gdis_sub,"\n");
    }
    
}

#' @export
popmap.out <- function(ofile,amp_euc,amp_bin,omit_popid){
    if(length(omit_popid)>0){
        ofi=sprintf("%s_euc%d_bin%d_omit%s_",ofile,amp_euc,amp_bin,paste0(omit_popid,collapse="-"));
    }else{
        ofi=sprintf("%s_euc%d_bin%d_",ofile,amp_euc,amp_bin);
    }

        for(k in 2:7){
            pngout(k,plotfile=paste0(ofi,k),outeps=T,density=700,geometry=700);
        }
    k=8;pngout(k,plotfile=paste0(ofi,k),outeps=T,density=700,geometry=1500);

}

#' @export
popmap.out.dev <- function(devid,ofile,amp_euc,amp_bin,omit_popid,density=700,geometry=700){
    if(length(omit_popid)>0){
        ofi=sprintf("%s_euc%d_bin%d_omit%s_",ofile,amp_euc,amp_bin,paste0(omit_popid,collapse="-"));
    }else{
        ofi=sprintf("%s_euc%d_bin%d_",ofile,amp_euc,amp_bin);
    }

            pngout(devid,plotfile=ofi,outeps=T,density=density,geometry=geometry);


}



##source("~/migseq/popmapper/popmapper0i_test4c_try_check.R");
##.ess.source('~/migseq/hirano_marimo_2021_08_04-2022_08_17/popmapper_all_marimo_utree_subdivide_nth_marimo.R', visibly = FALSE, output = TRUE)

max_depth <<- 0;

#' @export
subdata <- function(data,map,tagcid){
    cid=map$cid_sm;
    ms=(cid==tagcid);
    mss=(data$sample.info==tagcid);
    datab=list(dist2=data$dist2[ms,ms],
                   dist2p=data$dist2p[ms,ms],
                   ids2=data$ids2[ms],
                   popids2=data$popids2[ms],
                   popnames2=data$popnames2[ms],
                   sample_names2=data$sample_names2[ms],
                   sample_names1=data$sample_names1,
                   popid.legend=data$popid.legend,
                   popname.legend=data$popname.legend,
                   tab00=data$tab00,
                   tab1=data$tab1,
                   tab1p=data$tab1p,
                   ids1=data$ids1,
                   popids1=data$popids1,
                   runids1=data$runids1,
                   popids.si=data$popids.si,
                   popnames.si=data$popnames.si,
                   neu=data$neu,
                   neup=data$neup,
                   masks=data$masks,
                   mask=data$mask,
                   nrun=data$nrun,
                   runlabel=data$runlabel,
                   runid=data$runid,
                   sample.info=data$sample.info,
                   amp_euc=data$amp_euc,
                   amp_bin=data$amp_bin,
                   omit_popid=data$omit_popid,
                   sample_group_name=data$sample_group_name
                   );
        return(datab);
}




#' @export
popmap.find.sub <- function(data,map,param=param0,sample_group_name="gid"){
    if(flag_envs).ee.append("popmap.find.sub",environment());
    cid=map$cid_sm;
    ngroup=map$ngroup;
    maps.sub=NULL;
    for(i in 1:ngroup){
        maps.sub=c(maps.sub,list("none"));
        names(maps.sub)[i]=paste0("gid",i);
        
    }
    for(i in 1:ngroup){
        cat("\n###",i,"\n");
        
        cat("Sub-mapping for GroupID:",i,"\n");
        datab=subdata(data,map,i);
        sample_group_name1=paste0(sample_group_name,":",i);

        dev.set(6);
        plot_power_map(map$pcoa,map$powermap,map$powermap$cid_sm,ccer0=map$powermap$ccer0,ccer1=map$powermap$ccer1,labels=map$popids2,mask=i,main=sample_group_name1,param=param);

        flag_findsub=1;
        ##if(dim(datab$dist2)[1]<3){
        if(length(datab$dist2)< subdiv_edge_size^2){
            flag_findsub=0;
        }else{
            
            res= cmdscale(datab$dist2, k = 2,eig=TRUE);
            if(sum(res$eig>0)<2)flag_findsub=0;
        }
        if(flag_findsub==0){
            names(maps.sub)[i]=paste0(names(maps.sub)[i],".none");
        }
        
        if(flag_findsub==1){
        
            
     
        ##popmap.plotdata(map$pcoa,data,plot_mask=(map$cid_sm==i),col1="red",col2="black",tagid=sample_group_name1);
        
             
        maps.sub[[i]]=(popmap.find(datab,param=param,plot_data=0,sample_group_name=sample_group_name1));
        ##maps.sub[[i]]=c(maps.sub[[i]],list(data.sub=datab));
        maps.sub[[i]]=append(maps.sub[[i]],list(data.sub=datab));
        
        if(names(maps.sub[[i]])[1]=="pcoa"){
            names(maps.sub)[i]=paste0(names(maps.sub)[i],".none");
        }
    }
}

    
    return(maps.sub);
}



popmap.find.subs <- function(maps.sub,param=param0,sample_group_name="gid"){
    if(flag_envs).ee.append("popmap.find.subs",environment());
    maps.sub2=maps.sub;
    for(i in 1:length(maps.sub)){
        if(length(grep("none",names(maps.sub2)[i]))==0){
            
            maps.sub.sub=popmap.find.sub(maps.sub2[[i]]$data.sub,maps.sub2[[i]][[1]],param=param,sample_group_name=paste0(sample_group_name,":",i));
            ##map2=popmap.merge.sub(maps.sub2[[i]][[1]],maps.sub.sub);
            maps.sub2[[i]]=append(maps.sub2[[i]],list(maps.sub=maps.sub.sub));

            ##maps.sub2[[i]]=append(list(sub=map2),maps.sub2[[i]]);
            ##names(maps.sub2[[i]])[1]=paste0(names(maps.sub2[[i]])[2],".merge");
            
    }
    }
    return(maps.sub2);
}



popmap.find.subrec <- function(maps.sub,param=param0,depth=1,sample_group_name="gid"){
    if(flag_envs).ee.append("popmap.find.subrec",environment());    
    maps.sub2=popmap.find.subs(maps.sub,param=param,sample_group_name=sample_group_name);
    for(i in 1:length(maps.sub2)){
        not_none= (length(grep("none",names(maps.sub2)[i]))==0);
        not_merge= (length(grep("merge",names(maps.sub2)[i]))==0);
        if(not_none&&not_merge){
            cat("depth:",depth,"\n");
            if(max_depth < depth)max_depth <<- depth;
            print(paste("max depth:",max_depth));
            
            maps.sub3=popmap.find.subrec(maps.sub2[[i]]$maps.sub,param=param,depth=depth+1,sample_group_name=paste0(sample_group_name,":",i));
            maps.sub2[[i]]$maps.sub=maps.sub3;
        }
    }


 
        return(maps.sub2); 
}

#' @export
popmap.subdiv <- function(data,map,param=param0,noprint=FALSE){
    if(flag_envs).ee.append("popmap.subdiv",environment());
    if(class(map)=="popmaps"){
        maps.all=map;
        map=map[[1]];
    }else{
        maps.all=list(map,map$pcoa);
        names(maps.all)=c(paste0("ngroup",map$ngroup),"pcoa");
    }

    if(length(data$sample.info)==0)data$sample.info=map$sample.info;
    
    dev.set(2);
    ##plot.text(map);
    popmap.plot(map,peak=F,param=param);
    maps.sub=popmap.find.sub(data,map,param=param);
    maps.rec=popmap.find.subrec(maps.sub,param=param);
    ##maps.all=maps;
    
    maps.all$maps.sub=maps.rec;
    
    return(structure(maps.all,class="popmaps"));
}



########################################################################

popmap.merge.sub <- function(map,maps.sub0,depth=1,noprint=FALSE){
    if(flag_envs).ee.append("popmap.merge.sub",environment());
    ##map=maps[[1]];
    si=map$sample.info;

    si=data.frame(si,GroupID.sub=si$GroupID*(10^depth));
    cid_sm.sub=map$cid_sm*(10^depth);
    cid_sm=map$cid_sm;

    ngroup=map$ngroup;
    ##ngroup=max(cid_sm);
    
    for(i in 1:ngroup){
        if(length(grep("none",names(maps.sub0)[i]))==0){
            if(names(maps.sub0[[i]])[1]!="pcoa"){
                if(depth==1){
                    si.cid_sm.chi=(maps.sub0[[i]])[[1]]$sample.info$GroupID;
                    cid_sm.chi=(maps.sub0[[i]])[[1]]$cid_sm;
                }
                if(depth>1){
                    si.cid_sm.chi=(maps.sub0[[i]])[[1]]$sample.info$GroupID.sub;
                    cid_sm.chi=(maps.sub0[[i]])[[1]]$cid_sm.sub;
                    
                }

                cat("merge:")
                print(length(si$GroupID.sub));
                print(length(si.cid_sm.chi));
                si$GroupID.sub=si$GroupID.sub+si.cid_sm.chi;
                cid_sm.sub[cid_sm==i]=cid_sm.sub[cid_sm==i]+cid_sm.chi;
            }
        }
        }

    
    cidu=sort(unique(cid_sm.sub));
    cid=cid_sm.sub*0;
    for(i in 1:length(cidu)){
        cid[cid_sm.sub==cidu[i]]=i;        
    }

    si.cidu=sort(unique(si$GroupID.sub));
    si.cid=si$GroupID.sub;
    si.cid1=si$GroupID.sub*0;
    for(i in 1:length(si.cidu)){
       si.cid1[si.cid==si.cidu[i]]=i-as.integer(min(si.cidu)==0);
    }

    si=data.frame(si,GroupID.sub1=si.cid1);
    
    map$sample.info=si;
    map$cid_sm.sub=cid_sm.sub;

    map$cid_sm.merge=cid;
    ##cat("\n merge::",map$cid_sm.merge);
    map$ngroup=max(cid);
    if(!noprint){
        print(cid_sm.sub);
        print(cid);
    }
    return(map);
    
}




merge.rec <- function(maps.sub,param=param0,depth=3,noprint=FALSE){
    ##if(flag_envs).ee.append("find.subs",environment());
    maps.sub2=maps.sub;
    for(i in 1:length(maps.sub)){
        if(length(grep("none",names(maps.sub2)[i]))==0){
            ##maps.sub.sub=popmap.find.sub(maps.sub2[[i]]$data.sub,maps.sub2[[i]][[1]],param=param);
            maps.sub.sub=maps.sub[[i]]$maps.sub;
            if(depth>1){
                maps.sub.sub=merge.rec(maps.sub.sub,param=param,depth=depth-1,noprint=noprint);
                if(!noprint)print(paste("depth:",depth));
                map2=popmap.merge.sub(maps.sub2[[i]][[1]],maps.sub.sub,depth=depth-1,noprint=noprint);
                maps.sub2[[i]]=append(maps.sub2[[i]],list(maps.sub=maps.sub.sub));
                maps.sub2[[i]]=append(list(sub=map2),maps.sub2[[i]]);
                names(maps.sub2[[i]])[1]=paste0(names(maps.sub2[[i]])[2],".merge");
            }
    }
    }
    return(maps.sub2);
}


merge.map <- function(maps,maps.rec,depth=10,noprint=FALSE){
    if(flag_envs).ee.append("merge.map",environment());
    if(depth==0){
        map1=maps[[1]];
        map1$cid_sm.merge=map1$cid_sm;
        return(map1);
    }
    maps.rec1=merge.rec(maps.rec,depth=depth,noprint=noprint);
    map1=popmap.merge.sub(maps[[1]],maps.rec1,depth=depth,noprint=noprint);
    ##cat("\n ####:",map1$cid_sm.merge);
    if(!noprint){
        print("cid_sm.sub");
        map1$cid_sm.sub;
    }
    return(map1);
}

######################################################################

get.submaps <- function(maps.all,gid=NULL){
    if(length(gid)>0){
        com=sprintf("maps.t=maps.all$maps.sub$gid%d",gid[1]);
        if(length(gid)>1){
            for(i in 2:length(gid)){
                com=paste0(com,sprintf("$maps.sub$gid%d",gid[i]));
            }
        }
        
        eval(parse(text=com));

    }else{
        maps.t=maps.all;
        }
    return(maps.t);
}

#' @export
get.submap <- function(maps.all, gid=NULL,depth=NULL,noprint=FALSE,data=NULL){
    if(flag_envs).ee.append("get.submap",environment());
    
    if(length(gid)>0){
        com=sprintf("maps.t=maps.all$maps.sub$gid%d",gid[1]);
        if(length(gid)>1){
            for(i in 2:length(gid)){
                com=paste0(com,sprintf("$maps.sub$gid%d",gid[i]));
            }
        }
        
        eval(parse(text=com));
        dep=depth;
        if(length(dep)==0)dep=depth.maps(maps.t);


        map.t=merge.map(maps.t,maps.t$maps.sub,depth=dep,noprint=noprint);
        cat("depth for submap:", depth.maps(maps.t),"\n");
    }else{
        if(length(depth)==0)depth=depth.maps(maps.all);
        map.t=merge.map(maps.all,maps.all$maps.sub,depth=depth,noprint=noprint);
        cat("depth for maps:", depth.maps(maps.all),"\n");


    }

    if(length(data)>0){
        gdis.all=cmdist(data$dist2.all,map.t$cid_sm.merge,test=FALSE,flag_welch_t_test=param1$flag_welch_t_test,nperm=param1$nperm,flag_pforeach=0);
        colnames(gdis.all$gdis)=gid2sub(map.t);
        rownames(gdis.all$gdis)=gid2sub(map.t);
        map.t$gdist.all.sub=gdis.all;
    }

    
    return(map.t);

}

#' @export
plot.sub <- function(maps.all,gid=NULL,depth=NULL,labels="popid", main=NULL,noprint=FALSE,power=FALSE,peak=FALSE,main_add=TRUE,param=param0){
    if(length(depth)==0){
            depth_name="NULL";
    }else{
        depth_name=paste0(depth);
    }
    if(length(gid)==0){
            gid_name="all";
    }else{
        gid_name=paste0(gid,collapse=":");
    }
    
    if(length(main)==0)main=sprintf("gid: %s     depth: %s",gid_name,depth_name);
    if(class(labels)=="character")main=sprintf("%s    labels: %s\n",main,labels);
    
    if(power==FALSE){
        map.all=get.submap(maps.all,gid=gid,depth=depth,noprint=noprint);
        plot.text(map.all,labels=labels,main=main);
    }else{
        maps.t=get.submaps(maps.all,gid=gid);
        map.all=maps.t[[1]];
        popmap.plot(map.all,main=main,peak=peak,main_add=main_add,param=param);
    }
    
    invisible(map.all);
}


#####################################################################

depth.map <-function(map.t){
       
    cids=map.t$cid_sm.sub*10;
    if(length(cids)==0)return (0);
    dep=0;
    count=0;
    while(max(cids)>0){
        cids=cids%/%10;
        count=count+1;
        if(max(cids%%10)==0)dep=dep+1;
        ##cat(count,dep,"\n");
    }
    return (count-dep-1);
}
　

depth.maps <-function(maps.t){
 
    map.t=merge.map(maps.t,maps.t$maps.sub,depth=10,noprint=TRUE);
    return(depth.map(map.t));
}





gid2sub <- function(map.all){
    ma=sort(unique(map.all$cid_sm.merge));
    ma1=ma;
    if(length(map.all$cid_sm.sub)>0){
        ma1=ma*0;
        for(i in 1:length(ma1)){
            ma1[i]=map.all$cid_sm.sub[which(map.all$cid_sm.merge==ma[i])[1]];
        }
    }
    return(ma1);
}



plot.text <- function(map,labels=map$popids2,main="Subdivided-grouping map"){

    par(mar=c(5,5,5,6));
    par(xpd=T);
  if(length(labels)==1){
      if(labels=="gid.sub"){
          if(length(map$cid_sm.sub)>0){
              labels=map$cid_sm.sub;
          }else{
              print("No gid.sub found!");
              labels=map$popids2;
          }
      }else{
          lcase=c("gid","gid.merge","popid");
          labels0=list(gid=map$cid_sm,gid.merge=map$cid_sm.merge,popid=map$popids2);
          labels=labels0[[which(names(labels0)==labels)]];
      }
  }
    
    pals=rainbow(map$ngroup,v=1.0);
    ##pals=(rainbow(map$ngroup,end=0.7,v=0.8))
    xx0=map$pcoa$points[,1];
    yy0=map$pcoa$points[,2];
    xx=map$powermap$xx0;
    yy=map$powermap$yy0;
    zz=map$powermap$zz;

    dev.hold();
     plot(xx,yy,type="n",xlim=c(min(xx),max(xx)),ylim=c(min(yy),max(yy)),main=main,xaxs="i",yaxs="i",xlab="Axis 1", ylab="Axis 2");

    
    nlev=20;
    .filled.contour(xx,yy,zz,levels=seq(0,max(zz),,nlev),col=gray.colors(nlev,start=0.1,end=0.6));
    
    legend(par()$usr[2], par()$usr[4], legend = paste("group",seq(length(pals))), col = pals,pch=rep(15,length(pals)),cex=0.8,text.col="white",bg="#555555");

    text(x=xx0,y=yy0,col=pals[map$cid_sm.merge],labels=labels,cex=0.7,font=1);
    dev.flush();
}


plot.mask <- function(data,labels=data$popids2,mask=(data$ids2>0),col="red",cex=0.7){

    map=maps[[1]];tagid=2;mask=(map$cid_sm==tagid);
    labels=map$popids2[mask];
    xx0=maps$pcoa$points[mask,1];
    yy0=map$pcoa$points[mask,2];
    text(x=xx0,y=yy0,col="green",labels=labels,cex=0.7,font=1);
    
    
}

#' @export
plot_gdis <- function(map.all,gdis=map.all$gdist.all.sub$gdis,labels=NULL,main=NULL){
        cid=map.all$cid_sm.merge;
        
        re=cmdscale(gdis);
        popidd=map.all$popids2%/%10;
        pids=rep("",max(cid));
        npi=rep(0,max(cid));
        
        for(i in 1:max(cid)){
            pids[i]=paste(sort(unique(popidd[cid==i])),collapse=",");
            npi[i]=length(sort(unique(popidd[cid==i])));
        }
    ##plot(re,type="n",xlim=c(-10,25),ylim=c(-15,20));

        if(length(labels)==0)labels=pids;
        
    par(fg="black",bg="white");
    plot(re,type="n",xlab="Axis 1",ylab="Axis2",main=main);


    for(i in 2:ncol(gdis)){
        for(j in 1:(i-1)){
            lines(re[c(i,j),1],re[c(i,j),2],col="gray");
        }
    }

        
    text(re[,1],re[,2],labels=labels,col=rainbow(ncol(gdis),v=0.7),cex=1.5+0.5/npi);

        for(i in 2:ncol(gdis)){
        for(j in 1:(i-1)){
            text(mean(re[c(i,j),1]),mean(re[c(i,j),2]),labels=sprintf("%2.1f",gdis[i,j]),cex=0.8)
        }
        }
    }


##################################### for boot strap tree #####################

boot_dist <- function(tab){
    dist0=calc_dist(tab,amp_euc=1,amp_bin=1,flag_binary=2,flag_pforeach=1)$dist;
    ids1=data$ids1;       
    idsu=unique(ids1);    
    dist0u=sum_ind_wise(dist0,idsu);
    dist2=dist0u;
    ##ids2=idsu;
    ##colnames(dist2)=si$PopID[si$omit==0];
    ##rownames(dist2)=si$PopID[si$omit==0];
    colnames(dist2)=si$SampleName[si$omit==0];
    rownames(dist2)=si$SampleName[si$omit==0];
    
return(dist2);
}


pop_label <<- 0;
#' @export
boot_tree <- function(tab){
        dist0=calc_dist(tab,amp_euc=1,amp_bin=1,flag_binary=2,flag_pforeach=1)$dist;
    ids1=data$ids1;       
    idsu=unique(ids1);    
    dist0u=sum_ind_wise(dist0,idsu);
    dist2=dist0u;
        ##ids2=idsu;
        if(pop_label==1){
            colnames(dist2)=si$PopID[si$omit==0];
            rownames(dist2)=si$PopID[si$omit==0];
            }else{
                colnames(dist2)=data$sample_names2;
                rownames(dist2)=data$sample_names2;
            }
        tr=nj(dist2);
    
return(tr);
}

#' @export
calc_boot <- function(tab00,nsamp=100){
    nct=ncol(tab00);
##    tree=boot_tree(tab1);
    ##btree=list(NULL);
    for(i in 1:nsamp){
        print(i);
        ##lis=order(rnorm(ncol(tab00)));
        ##tabb=tab00[,lis[1:(length(lis)%/%ndiv)]];
        
        lis=ceiling(runif(nct)*nct);
        tabb=tab00[,lis];
        tr1p=boot_tree(tabb);
        if(i==1)btree=list(tr1p);
        if(i>1)btree=c(btree,list(tr1p));
    }
        return(btree);
}



#' @export
plot.btree1 <- function(btree,map.all,ccid=map.all$cid_sm,p_boot=0.95,power=0.5,type="u",cex_boot=0.5,cex=1.0,ecol="white",fg="white",bg="black",outgroup=NULL,sample_names=NULL){
    ##if(flag_envs).ee.append("plot.btree1",environment());    

   
    btree1=btree;

    ctree=consensus(btree1,p=p_boot);
##    if(length(outgroup)>0)ctree=RootTree(ctree,outgroup);
    ctree1=ctree;

    
    ctree1$tip.label=map.all$popids2;
    con=consensus.edges(btree1, consensus.tree=ctree,method="least.squares");
    ##con=consensus.edges(btree, consensus.tree=ctree,method="mean.edge");
    ctree1$edge.length=con$edge.length;

    boot = prop.clades(ctree,btree1);
    
    ##ccid=map.all$cid_sm;
    if(length(map.all$pals.merge)>0){
        pals=map.all$pals.merge;
    }else{
        pals=(rainbow(max(ccid),end=0.7,v=0.8))
    }
    
    ecols=rep(ecol,length(con$edge.length));
    ecols[con$edge.length<=0]="blue";
    lmin=min(con$edge.length[con$edge.length>0]);
    ct2=ctree1;
    ct2$edge.length[ct2$edge.length<=0]=lmin;
    ct2$edge.length=(ct2$edge.length)^power;


    if(type=="p"){
        ##if(length(sample_names)==0)sample_names=data$sample_names2;
        ##ct2$tip.label=paste0(ct2$tip.label," -- ",sample_names);
        if(length(sample_names)>0)ct2$tip.label=sample_names;
        }
    
    ##par(bg="#999999",fg="black");
par(bg=bg,fg=fg);


if(length(outgroup)>0){
    print(outgroup);
    ogid=rep(0,length(outgroup));
    for(i in 1:length(outgroup))ogid[i]=as.integer(grep(outgroup[i],ct2$tip.label));
    ogid=as.integer(ogid);
    
    ## boot = prop.clades(ct2,btree1);
    ct2_root=root(ct2,ogid);
    boot_root = prop.clades(root(ctree,ogid),btree1);
    
    plot(ct2_root,direction="u",tip.color=pals[ccid],type=type,edge.color=ecols,cex=cex);

    
    drawSupportOnEdges(boot_root,cex=cex_boot);
}else{
    
    plot(ct2,direction="u",tip.color=pals[ccid],type=type,edge.color=ecols,cex=cex);

    drawSupportOnEdges(boot,cex=cex_boot);
    }
}
