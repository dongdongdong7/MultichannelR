library(tidyverse)
# mix3

# ****** 0. Initialize ****** #
{
  # 0.1 Configuration parameter
  {
    deltaMz <- 2.0125 # Da
    deltaRt <- 30 # s
    tolMz1 <- 0.01 # Da
    tolMz2 <- 0.01 # Da
    
    deltaIsotope1 <- 1.0033
    deltaIsotope2 <- deltaIsotope1 / 2
    deltaIsotope3 <- deltaIsotope1 / 3
    
    channelNumber <- 6
    Q3_A1 <- c(171.1038, 172.1101, 173.1165, 174.1229, 175.1291, 176.1351)
    Q3_A2 <- c(198.1276, 200.1398, 202.1527, 204.1652, 206.1776, 208.1899)
    Q3_A3 <- c(381.1267, 383.1393, 385.1581, 387.1644, 389.1770, 391.1895)
    Q3_P1 <- c(199.1351, 201.1482, 203.1604, 205.1734, 207.1857, 209.1986)
    Q3 <- c(Q3_A1, Q3_A2, Q3_A3, Q3_P1)
    
    cwp_centroid <- xcms::CentWaveParam(snthresh = 0, noise = 0, ppm = 10,
                                        peakwidth = c(1, 10), firstBaselineCheck = FALSE)
  }
  
  # 0.2 Load data
  {
    res_dir <- "D:/fudan/Projects/2024/MultichannelR/Progress/build_package/MultichannelR/"
    dir_path <- "D:/fudan/Projects/2024/MultichannelR/Progress/build_package/MultichannelR/test_data/mix6/"
    files_path <- list.files(dir_path, pattern = ".mzML")
    files_path <- paste0(dir_path, files_path)
    
    MSnbaseData <- MSnbase::readMSData(files_path, mode = "onDisk")
    sps <- Spectra::Spectra(files_path, backend = Spectra::MsBackendMzR())
    
    sampleNames <- MSnbase::pData(MSnbaseData)$sampleNames
    sampleNumber <- length(sampleNames)
  }
  
  # 0.3 Load Table
  {
    # dataTable
  }
  
  # 0.4 Function
  {
    # 1. 手动获取一个peak的2级质谱
    # 1.1 手动获得一个peak的所有2级质谱
    # sps 是一个Spectra对象,包含要获取的2级质谱信息
    # peaksInfo 是由as.data.frame(xcms::chromPeaks(obj))获取而来
    # peakName 是peak的名字,peaksInfo的行名 # 这个deltaMz应该为tolMz
    GetSpectra2 <- function(sps, peaksInfo, peakName, tolMz = 0.001, deltaRt = 0){
      peak <- peaksInfo[peaksInfo$peak_id == peakName, ]
      sp2 <- sps %>% 
        Spectra::filterMsLevel(2) %>% 
        Spectra::filterPrecursorMzRange(c(peak$mzmin - tolMz, peak$mzmax + tolMz)) %>% 
        Spectra::filterRt(c(peak$rtmin - deltaRt, peak$rtmax + deltaRt))
      
      print(paste0("MS2 number: ", length(sp2)))
      
      if(length(sp2) == 0){
        sp2 <- NULL
      }else{
        sp2$peak_id <- peakName
      }
      
      return(sp2)
    }
    # 1.2 一个peak的多个二级质谱中离peak RT最近的那一个
    GetSpectra2_closeRT <- function(sps, peaksInfo, peakName, tolMz = 0.001, deltaRt = 0){
      peak <- peaksInfo[peaksInfo$peak_id == peakName, ]
      rt <- peak$rt
      sp2 <- GetSpectra2(sps, peaksInfo, peakName, tolMz, deltaRt)
      if(length(sp2) > 1){
        idx <- which.min(abs(Spectra::rtime(sp2) - rt))
        sp2 <- sp2[idx]
      }
      return(sp2)
    }
    # 1.3 将一个质谱标准化
    # sp2是一个单独的spectra
    standardizeSpectra <- function(sp2){
      
      if(is.null(sp2)){
        return(sp2)
      }
      
      norm_fun <- function(z, ...) {
        z[, "intensity"] <- z[, "intensity"] /
          max(z[, "intensity"], na.rm = TRUE) * 100
        z
      }
      
      sp2 <- Spectra::filterIntensity(sp2, intensity = max(Spectra::intensity(sp2)) * 0.05)
      sp2 <- Spectra::addProcessing(sp2, FUN = norm_fun)
      
      return(sp2)
    }
    # 2. 绘制质谱
    # 2.1 绘制一张spectra
    # sp2 就是Spectra package中的一个spectra
    # num 就是标记int为前num的singnal
    plot_spectra <- function(sp2, num = 5){
      mz <- unlist(Spectra::mz(sp2))
      int <- unlist(Spectra::intensity(sp2))
      df <- tibble(mz = mz, int = int) %>% 
        arrange(desc(int))
      n <- nrow(df)
      if(num > n){
        num <- n
      }
      label_point <- c(df$int[1:num], rep(NA, n - num))
      label_text <- as.character(c(sprintf("%.4f", df$mz[1:num]), rep(NA, n-num)))
      df$label_point <- label_point
      df$label_text <- label_text
      
      p <- ggplot(df, aes(x = mz, y = int)) +
        geom_segment(aes(x = mz, xend = mz, y = 0, yend = int), color = "grey") + 
        geom_point(aes(x = mz, y = label_point), color = "orange", size = 4) + 
        ggrepel::geom_text_repel(aes(label = label_text), size = 4, vjust = -1, min.segment.length = Inf) + 
        labs(title = sp2$peak_id, subtitle = paste0("Retention time: ", sprintf("%.4f", Spectra::rtime(sp2)))) + 
        theme_classic()
      return(p)
    }
    # 3. 获取特征碎片和分类
    # sp2 为2级质谱
    # tolMz 为在这个tolMz之中的mz被认为是同一离子
    # num为搜索前num个离子
    GetSpecialFragment <- function(sp2, tolMz = 0.01, num = 10){
      Q1 <- Spectra::precursorMz(sp2)
      mz <- unlist(Spectra::mz(sp2))
      int <- unlist(Spectra::intensity(sp2))
      df <- tibble(mz = mz, int = int) %>% 
        arrange(desc(int)) %>% 
        filter(!near(mz, Q1, tol = tolMz))
      n <- nrow(df)
      if(n == 0){
        return(list(NA, NA))
      }
      if(num > n){
        num <- n
      }
      for(i in 1:nrow(df)){
        Q3 <- df$mz[i]
        logical_A1 <- near(Q3, Q3_A1, tol = tolMz)
        logical_A2 <- near(Q3, Q3_A2, tol = tolMz)
        logical_A3 <- near(Q3, Q3_A3, tol = tolMz)
        logical_P1 <- near(Q3, Q3_P1, tol = tolMz)
        if(any(logical_A1)){
          Q3 <- Q3_A1[which(logical_A1 == TRUE)]
          class <- "A1"
          break
        }else if(any(logical_A2)){
          Q3 <- Q3_A2[which(logical_A2 == TRUE)]
          class <- "A2"
          break
        }else if(any(logical_P1)){
          Q3 <- Q3_P1[which(logical_P1 == TRUE)]
          class <- "P1"
          break
        }else if(any(logical_A3)){
          Q3 <- Q3_A3[which(logical_A3 == TRUE)]
          class <- "A3"
          break
        }else{
          Q3 <- NA
          class <- NA
          next
        }
      }
      return(list(Q3, class))
    }
    # 4. 寻找同位素
    # 4.1 输入是一个peakInfo
    # peakInfo <- peakGroupList[[236]][1, ]
    # peaksInfo_ALL 是未经过筛选的全部peaksInfo
    findIso <- function(peakInfo, peaksInfo_ALL, tolMz = 0.001, tolRt = 1){
      
      deltaIsotope1 <- 1.0033
      deltaIsotope2 <- deltaIsotope1 / 2
      deltaIsotope3 <- deltaIsotope1 / 3
      deltaIsotope4 <- deltaIsotope1 / 4
      
      MZ <- peakInfo$mz
      RT <- peakInfo$rt
      mz_isotope1 <- MZ + deltaIsotope1
      mz_isotope2 <- MZ + deltaIsotope2
      mz_isotope3 <- MZ + deltaIsotope3
      mz_isotope4 <- MZ + deltaIsotope4
      
      #browser()
      tmp1 <- peaksInfo_ALL %>% 
        filter(near(mz, mz_isotope1, tolMz) & near(rt, RT, tolRt))
      idx1 <- which(tmp1$maxo < peakInfo$maxo)
      if(length(idx1) == 0){
        tmp1 <- tmp1[0, ]
      }else{
        tmp1 <- tmp1[idx1, ]
      }
      
      tmp2 <- peaksInfo_ALL %>% 
        filter(near(mz, mz_isotope2, tolMz) & near(rt, RT, tolRt))
      idx2 <- which(tmp2$maxo < peakInfo$maxo)
      if(length(idx2) == 0){
        tmp2 <- tmp2[0, ]
      }else{
        tmp2 <- tmp2[idx2, ]
      }
      
      tmp3 <- peaksInfo_ALL %>% 
        filter(near(mz, mz_isotope3, tolMz) & near(rt, RT, tolRt))
      idx3 <- which(tmp3$maxo < peakInfo$maxo)
      if(length(idx3) == 0){
        tmp3 <- tmp3[0, ]
      }else{
        tmp3 <- tmp3[idx3, ]
      }
      
      tmp4 <- peaksInfo_ALL %>% 
        filter(near(mz, mz_isotope4, tolMz) & near(rt, RT, tolRt))
      idx4 <- which(tmp4$maxo < peakInfo$maxo)
      if(length(idx4) == 0){
        tmp4 <- tmp4[0, ]
      }else{
        tmp4 <- tmp4[idx4, ]
      }
      
      if(nrow(tmp4) >= 1){
        return(4)
      }else if(nrow(tmp3) >= 1){
        return(3)
      }else if(nrow(tmp2) >= 1){
        return(2)
      }else if(nrow(tmp1) >= 1){
        return(1)
      }else{
        return(NA)
      }
    }
    # 5. 输入specialFragment的数值, 返回他的sampleIdx
    GetSampleIdx <- function(specialFragment){
      Q3_A1 <- c(171.1038, 172.1101, 173.1165, 174.1229, 175.1291, 176.1351)
      Q3_A2 <- c(198.1276, 200.1398, 202.1527, 204.1652, 206.1776, 208.1899)
      Q3_A3 <- c(381.1267, 383.1393, 385.1581, 387.1644, 389.1770, 391.1895)
      Q3_P1 <- c(199.1351, 201.1482, 203.1604, 205.1734, 207.1857, 209.1986)
      Q3 <- c(Q3_A1, Q3_A2, Q3_A3, Q3_P1)
      
      if(length(specialFragment) != 1){
        stop("SpecialFragment's length is not 1")
      }
      if(!specialFragment %in% Q3){
        return(NA)
      }else{
        logical_A1 <- specialFragment %in% Q3_A1
        logical_A2 <- specialFragment %in% Q3_A2
        logical_A3 <- specialFragment %in% Q3_A3
        logical_P1 <- specialFragment %in% Q3_P1
        
        if(logical_A1){
          sampleIdx <- which(Q3_A1 %in% specialFragment)
        }else if(logical_A2){
          sampleIdx <- which(Q3_A2 %in% specialFragment)
        }else if(logical_A3){
          sampleIdx <- which(Q3_A3 %in% specialFragment)
        }else if(logical_P1){
          sampleIdx <- which(Q3_P1 %in% specialFragment) 
        }else{
          sampleIdx <- NA
        }
      }
      return(sampleIdx)
    }
    # 6. 输入一个peakGroup, 返回peakGroupTibble
    # peakGroup <- peakGroupList_list[[1]][[309]]
    GetpeakGroupTibble <- function(peakGroup){
      peaksNum <- nrow(peakGroup)
      tagNum <- peakGroup$tagNum
      length1 <- length(which(tagNum == 1))
      length2 <- length(which(tagNum == 2))
      length3 <- length(which(tagNum == 3))
      length4 <- length(which(tagNum == 4))
      lengthNA <- length(which(is.na(tagNum)))
      if(length1 > length2 & length1 > length3 & length1 > length4 & length1 != 0){
        tagNum <- 1
      }else if(length2 > length1 & length2 > length3 & length2 > length4 & length2 != 0){
        tagNum <- 2
      }else if(length3 > length1 & length3 > length2 & length3 > length4 & length3 != 0){
        tagNum <- 3
      }else if(length4 > length1 & length4 > length2 & length4 > length3 & length4 != 0){
        tagNum <- 4
      }else{
        tagNum <- NA
      }
      
      peak <- peakGroup[1, ]
      peakMz <- peak$mz
      peakRt <- peak$rt
      sample <- peak$sample
      class <- peak$class[[1]]
      peakChannelIdx <- peak$channelIdx[[1]]
      
      reagent <- MetaboCoreUtils::standardizeFormula("C14H15NO2S")
      reagent_mz <- MetaboCoreUtils::formula2mz(reagent, adduct = "[M+H]+")[[1]]
      reagent_mass <- MetaboCoreUtils::mz2mass(reagent_mz)[[1]]
      
      reagent_mz <- reagent_mz + (peakChannelIdx - 1) * 2
      
      # 如果知道了tagNum的数量, 我们可以计算衍生化前的mass
      if(!is.na(tagNum)){
        if(tagNum == 1){
          mass <- peakMz - reagent_mz
        }else if(tagNum == 2){
          mass <- peakMz * 2 - reagent_mz * 2
        }else if(tagNum == 3){
          mass <- peakMz * 3 - reagent_mz * 3
        }else if(tagNum == 4){
          mass <- peakMz * 4 - reagent_mz * 4
        }
      }else{
        mass <- NA
      }
      
      peakGroupTibble <- tibble(mass = mass, rt = peakRt, peaksNum = peaksNum, sample = sample, tagNum = tagNum, class = class, peakGroup = list(peakGroup))
      return(peakGroupTibble)
    }
    # 7. fill PeaksGroupTibbles
    # peaksTibble <- ida_data_centroid_peaksTibble
    # sps <- sps_centroid
    FillPeaksGroup <- function(PeaksGroupTibbles, sps, peaksTibble, tolRt_fill = 10){
      tf <- PeaksGroupTibbles %>% 
        filter(peaksNum != 6 & peaksNum >= 3)
      for(i in 1:nrow(tf)){
        PeaksGroupTibble <- tf[i, ]
        peaksNum <- PeaksGroupTibble$peaksNum
        PeaksGroup <- PeaksGroupTibble$peaksGroup[[1]]
        sampleIdxNew <- which(!c(1, 2, 3, 4, 5, 6) %in% PeaksGroup$sampleIdx)
        standardPeak <- PeaksGroup[1, ]
        standardIdx <- standardPeak$sampleIdx
        rt_mean <- mean(PeaksGroup$rt)
        for(j in 1:length(sampleIdxNew)){
          idx <- sampleIdxNew[j]
          if(idx > standardIdx){
            mz_new <- standardPeak$mz + deltaMz * (idx - standardIdx)
          }else{
            mz_new <- standardPeak$mz - deltaMz * (standardIdx - idx)
          }
          new_peak <- peaksTibble %>% 
            filter(near(mz, mz_new, tolMz1) & near(rt, rt_mean, tolRt_fill))
          if(nrow(new_peak) == 0){
            next
          }else if(nrow(new_peak) > 1){
            new_peak <- new_peak %>% 
              filter(maxo == max(maxo))
          }
          sp2 <- new_peak$spectra[[1]]
          if(!is.null(sp2)){
            new_peak$precursorMz <- precursorMz(sp2)
            tmp <- GetSpecialFragment(sp2)
            new_peak$specialFragment <- tmp[[1]][1]
            new_peak$class <- tmp[[1]][2]
            new_peak$sampleIdx <- idx
            new_peak$tagNum <- NA
          }else{
            new_peak$precursorMz <- NA
            new_peak$specialFragment <- NA
            new_peak$class <- NA
            new_peak$sampleIdx <- idx
            new_peak$tagNum <- NA
          }
          PeaksGroup <- rbind(PeaksGroup, new_peak)
        }
        PeaksGroup <- PeaksGroup %>% arrange(mz)
        tf[i, ] <- GetPeaksGroupTibble(PeaksGroup)
      }
      return(tf)
    }
    # 8. Get information about mz, rt and class of peak group on channel1
    # peakgroup <- peakGroupTibble_list[[1]][42, ]$peakGroup[[1]]
    getInfo_channel1 <- function(peakgroup, deltaMz = 2.0125){
      mz <- as.vector(peakgroup$mz[1] - (peakgroup$channelIdx[1] - 1) * deltaMz)
      rt <- as.vector(peakgroup$rt[1])
      class <- as.vector(peakgroup$class[1])
      df <- as.data.frame(matrix(NA, ncol = 3))
      colnames(df) <- c("mz", "rt", "class")
      df$mz <- mz;df$rt <- rt;df$class <- class
      return(df)
    }
  }
}
# ****** 0. Initialize ****** #

# ****** 1. Feature picking ****** #
{
  # 1.1 CentWave pick peaks
  {
    start_time <- Sys.time()
    MSnbaseDataPeaks <- xcms::findChromPeaks(MSnbaseData, param = cwp_centroid)
    end_time <- Sys.time()
    print(end_time - start_time)
    MSnbaseDataPeaksInfo <- as_tibble(xcms::chromPeaks(MSnbaseDataPeaks), rownames = "peak_id")
    peaksInfo_ALL <- MSnbaseDataPeaksInfo # 这个是最原始的peaksInfo, 不经过过滤
    # 进行过滤
    MSnbaseDataPeaksInfo <- MSnbaseDataPeaksInfo %>% 
      filter(maxo >= 1000)
    nrow(MSnbaseDataPeaksInfo) # 11010
  }
  
  # 1.2 Classify different samples
  {
    MSnbaseDataPeaksInfoList <- list()
    for(i in 1:sampleNumber){
      MSnbaseDataPeaksInfoList[[i]] <- MSnbaseDataPeaksInfo %>% filter(sample == i)
    }
    names(MSnbaseDataPeaksInfoList) <- sampleNames
    spsList <- list()
    for(i in 1:sampleNumber){
      spsList[[i]] <- sps %>% Spectra::filterDataStorage(unique(Spectra::dataStorage(sps))[i])
    }
    names(spsList) <- sampleNames
    peaksInfo_ALLList <- list()
    for(i in 1:sampleNumber){
      peaksInfo_ALLList[[i]] <- peaksInfo_ALL %>% filter(sample == i)
    }
    names(peaksInfo_ALLList) <- sampleNames
  }
  
  # 1.3 MS2 is assigned to each peak
  {
    peakSpectraList_list <- list()
    for(i in 1:sampleNumber){
      peaksInfo <- MSnbaseDataPeaksInfoList[[i]]
      peakSpectraList <- lapply(peaksInfo$peak_id, function(x){
        sp2 <- GetSpectra2_closeRT(spsList[[i]], peaksInfo, x)
        sp2 <- standardizeSpectra(sp2)
        return(sp2)
      })
      names(peakSpectraList) <- peaksInfo$peak_id
      peakSpectraList_list[[i]] <- peakSpectraList
    }
    names(peakSpectraList_list) <- sampleNames
    for(i in 1:sampleNumber){
      MSnbaseDataPeaksInfoList[[i]]$spectra <- peakSpectraList_list[[i]]
      logical_NULL <- map_vec(MSnbaseDataPeaksInfoList[[i]]$spectra, function(x) {!is.null(x)})
      MSnbaseDataPeaksInfoList[[i]] <- MSnbaseDataPeaksInfoList[[i]] %>% 
        filter(logical_NULL) %>% 
        arrange(mz)
    }
  }
}
# ****** 1. Feature picking ****** #

# ****** 2. Feature Grouping ****** #
{
  # 2.1 Each MS2 is assigned a parent ion
  {
    for(i in 1:sampleNumber){
      precursorMz_vec <- map_vec(MSnbaseDataPeaksInfoList[[i]]$spectra, function(x) {Spectra::precursorMz(x)})
      MSnbaseDataPeaksInfoList[[i]]$precursorMz <- precursorMz_vec
    }
  }
  
  # 2.2 Assign a feature fragment to each MS2 and state its class
  {
    for(i in 1:sampleNumber){
      tmpList <- lapply(MSnbaseDataPeaksInfoList[[i]]$spectra, function(x) {GetSpecialFragment(x)})
      specialFragment_vec <- map_vec(tmpList, function(x) {pluck(x, 1)})
      class_vec <- map_vec(tmpList, function(x) {pluck(x, 2)})
      MSnbaseDataPeaksInfoList[[i]]$specialFragment <- specialFragment_vec
      MSnbaseDataPeaksInfoList[[i]]$class <- class_vec
      MSnbaseDataPeaksInfoList[[i]] <- MSnbaseDataPeaksInfoList[[i]] %>% 
        filter(!is.na(specialFragment))
    }
  }
  
  # 2.3 Allocate channelIdx for each peak
  {
    for(i in 1:sampleNumber){
      channelIdx_vec <- map_vec(MSnbaseDataPeaksInfoList[[i]]$specialFragment, function(x) {GetSampleIdx(x)})
      MSnbaseDataPeaksInfoList[[i]]$channelIdx <- channelIdx_vec
    }
  }
  
  # 2.4 Find the isotope peak for each peak and measure the tag count
  {
    for(i in 1:sampleNumber){
      tagNum_vec <- map_vec(1:nrow(MSnbaseDataPeaksInfoList[[i]]), function(x) findIso(MSnbaseDataPeaksInfoList[[i]][x, ], peaksInfo_ALLList[[i]]))
      MSnbaseDataPeaksInfoList[[i]]$tagNum <- tagNum_vec
    }
  }
  
  # 2.5 Peak grouping
  {
    peakGroupList_list <- list()
    for(i in 1:sampleNumber){
      peakGroupList <- list()
      peaksTibble2 <- MSnbaseDataPeaksInfoList[[i]]
      peaksTibble3 <- peaksTibble2
      while(nrow(peaksTibble3) != 0){
        peaksTmp <- peaksTibble3[1, ]
        peak <- peaksTibble3[1, ]
        cpid <- peak$peak_id
        mz <- peak$mz
        rt <- peak$rt
        class <- peak$class[[1]]
        q3 <- peak$specialFragment[[1]]
        delta_mz <- peaksTibble3$mz - mz
        delta_rt <- abs(peaksTibble3$rt - rt)
        delta_q3 <- peaksTibble3$specialFragment - q3
        class_vec <- peaksTibble3$class
        start_loop <- peak$channelIdx
        delete_idx <- c(1)
        if(class == "A1"){
          deltaQ3 <- deltaMz / 2
        }else{
          deltaQ3 <- deltaMz
        }
        
        if(start_loop == 6){
          # do nothing
        }else{
          end_loop <- 6 - start_loop
          for(j in 1:end_loop){
            idx <- which(near(delta_mz, deltaMz * j, tol = tolMz1) & 
                           near(delta_rt, 0, tol = deltaRt) & 
                           class_vec %in% class & 
                           near(delta_q3, deltaQ3 * j, tol = tolMz2))
            if(length(idx) == 0){
              next
            }else if(length(idx) > 1){
              idx <- idx[which.max(peaksTibble3[idx, ]$intb)]
            }
            peaksTmp <- rbind(peaksTmp, peaksTibble3[idx, ])
            delete_idx <- c(delete_idx, idx)
          }
        }
        peaksTibble3 <- peaksTibble3[-delete_idx, ]
        peakGroupList <- append(peakGroupList, list(peaksTmp))
      }
      peakGroupList_list[[i]] <- peakGroupList
    }
    names(peakGroupList_list) <- sampleNames 
  }
  
  # 2.6 peakGroupList to peakGroupTibble
  {
    peakGroupTibble_list <- list()
    for(i in 1:sampleNumber){
      peakGroupTibble <- lapply(peakGroupList_list[[i]], function(x) GetpeakGroupTibble(x))
      peakGroupTibble <- list_rbind(peakGroupTibble)
      peakGroupTibble <- peakGroupTibble %>% 
        filter(!is.na(tagNum)) %>% 
        filter(mass > 0) %>% 
        filter(peaksNum >=3) %>% 
        arrange(mass)
      peakGroupTibble_list[[i]] <- peakGroupTibble
    }
    names(peakGroupTibble_list) <- sampleNames
  }
  
  # 2.7 rearrange iso
  {
    peakGroupTibble_list_back <- peakGroupTibble_list
    peakGroupTibble_list <- peakGroupTibble_list_back
    for(i in 1:length(peakGroupTibble_list)){
      peakGroupTibble <- peakGroupTibble_list[[i]]
      tmp <- lapply(1:nrow(peakGroupTibble), function(x){
        getInfo_channel1(peakGroupTibble[x, ]$peakGroup[[1]])
      })
      tmp <- list_rbind(tmp)
      isoIdx <- c()
      loop_idx <- 1:nrow(tmp)
      for(j in loop_idx){
        if(j %in% isoIdx) next
        mz1 <- tmp[j, ]$mz
        rt1 <- tmp[j, ]$rt
        class1 <- tmp[j, ]$class
        isoIdx_tmp <- which(near(mz1, tmp$mz, tol = 0.02) & near(rt1, tmp$rt, tol = 25) & class1 == tmp$class)
        #if(any(isoIdx_tmp %in% isoIdx)) next
        if(length(isoIdx_tmp) > 1){
          isoIdx <- c(isoIdx, isoIdx_tmp)
          peakGroupList_iso <- peakGroupTibble[isoIdx_tmp, ]$peakGroup
          for(k in 1:channelNumber){
            tmpList <- lapply(1:length(peakGroupList_iso), function(x){
              tb <- peakGroupList_iso[[x]][which(peakGroupList_iso[[x]]$channelIdx == k), ]
              if(nrow(tb) == 0) tb[1, ] <- NA
              return(tb)
            })
            tmpTibble <- list_rbind(tmpList)
            tmpTibble <- tmpTibble %>% arrange(rt)
            for(l in 1:length(isoIdx_tmp)){
              if(nrow(peakGroupList_iso[[l]][which(peakGroupList_iso[[l]]$channelIdx == k), ]) == 1){
                peakGroupList_iso[[l]][which(peakGroupList_iso[[l]]$channelIdx == k), ] <- tmpTibble[l, ]
              }else if(nrow(peakGroupList_iso[[l]][which(peakGroupList_iso[[l]]$channelIdx == k), ]) == 0){
                peakGroupList_iso[[l]] <- rbind(peakGroupList_iso[[l]], tmpTibble[l, ])
              }
              peakGroupList_iso[[l]] <- peakGroupList_iso[[l]] %>% 
                filter(!is.na(peak_id)) %>% 
                arrange(channelIdx)
            }
          }
          tmpList <- lapply(peakGroupList_iso, function(x) GetpeakGroupTibble(x))
          tmpTb <- list_rbind(tmpList)
          peakGroupTibble[isoIdx_tmp, ] <- tmpTb
        }
      }
      peakGroupTibble_list[[i]] <- peakGroupTibble
    } 
  }
}
# ****** 2. Feature Grouping ****** #

# ****** 3. Feature aligning ****** #
{
  # 3.1 Align
  {
    ref_idx <- 1
    all_idx <- 1:sampleNumber
    other_idx <- setdiff(all_idx, ref_idx)
    ref_peakGroupTibble <- peakGroupTibble_list[[ref_idx]]
    alignedGroupList <- list()
    for(i in 1:nrow(ref_peakGroupTibble)){
      refGroup <- ref_peakGroupTibble[i, ]
      alignedGroup <- refGroup
      ref_mass <- refGroup$mass
      ref_rt <- refGroup$rt
      ref_tagNum <- refGroup$tagNum
      ref_class <- refGroup$class
      for(j in other_idx){
        tmp_peakGroupTibble <- peakGroupTibble_list[[j]]
        tmp_mass <- tmp_peakGroupTibble$mass
        tmp_rt <- tmp_peakGroupTibble$rt
        tmp_tagNum <- tmp_peakGroupTibble$tagNum
        tmp_class <- tmp_peakGroupTibble$class
        
        idx <- which(near(ref_mass, tmp_mass, tol = tolMz1) &
                       near(ref_rt, tmp_rt, tol = 5) &
                       ref_tagNum == tmp_tagNum &
                       ref_class == tmp_class)
        if(length(idx) >= 1){
          alignedGroup <- rbind(alignedGroup, tmp_peakGroupTibble[idx, ])
        }
      }
      if(nrow(alignedGroup) != 1){
        alignedGroupList <- append(alignedGroupList, list(alignedGroup))
      }
    }
  }
  
  # 3.2 The number is not equal to the number of samples being filtered
  {
    # logical_delete <- map_vec(alignedGroupList, function(x){
    #   if(nrow(x) == sampleNumber){
    #     return(TRUE)
    #   }else{
    #     return(FALSE)
    #   }
    # })
    # alignedGroupList <- alignedGroupList[logical_delete]
  }
}
# ****** 3. Feature aligning ****** #

# ****** 4. Metabolite Identification ****** #
{
  # 4.1 Load database
  {
    AmidoLibrary <- read_tsv("D:/fudan/Projects/2024/MultichannelR/Progress/build_package/MultichannelR/Library/AmidoLibrary.txt")
    AmidoLibrary <- AmidoLibrary %>% 
      filter(Source == "HMDB")
    PheHydroLibrary <- read_tsv("D:/fudan/Projects/2024/MultichannelR/Progress/build_package/MultichannelR/Library/PheHydroLibrary.txt")
    PheHydroLibrary <- PheHydroLibrary %>% 
      filter(Source == "HMDB")
  }
  
  # 4.2 Mass match
  {
    resList <- list()
    for(i in 1:length(alignedGroupList)){
      feature_mass <- mean(alignedGroupList[[i]]$mass)
      feature_class <- alignedGroupList[[i]]$class[1]
      if(feature_class %in% c("A1", "A2", "A3")){
        res <- AmidoLibrary %>% 
          filter(near(feature_mass, Mass, tol = 0.05))
      }else if(feature_class %in% c("P1")){
        res <- PheHydroLibrary %>% 
          filter(near(feature_mass, Mass, tol = 0.05))
      }
      resList <- append(resList, list(res))
    }
  }
}
# ****** 4. Metabolite Identification ****** #

# ****** 5. Result Report ****** #
{
  # 5.1 创建report表格
  {
    report <- as.data.frame(matrix(NA, ncol = 25, nrow = 0))
    sampleIdx_channelIdx <- unlist(map(1:sampleNumber, function(x){
      paste0(x, "_", 1:channelNumber)
    }))
    colnames(report) <- c("Mass", "Rt", "Formula", "Class", "TagNumber", "Metabolites", "ID",
                          sampleIdx_channelIdx)
    report <- as_tibble(report)
    peakGroupNumber <- length(alignedGroupList)
    for(i in 1:peakGroupNumber){
      report_tmp <- report[0, ]
      report_tmp[1, ] <- NA
      alignedGroup <- alignedGroupList[[i]]
      res <- resList[[i]]
      Mass <- mean(alignedGroup$mass)
      Rt <- mean(alignedGroup$rt)
      Class <- unique(alignedGroup$class)[1]
      TagNumber <- unique(alignedGroup$tagNum)[1]
      if(nrow(res) == 0){
        Formula <- NA
        Metabolites <- NA
        ID <- NA
      }else{
        Formula <- paste0(unique(res$Formula), collapse = "|")
        Metabolites <- paste0(res$Name, collapse = "|")
        ID <- paste0(res$ID, collapse = "|")
      }
      report_tmp[1, "Mass"] <- Mass
      report_tmp[1, "Rt"] <- Rt
      report_tmp[1, "Class"] <- Class
      report_tmp[1, "TagNumber"] <- TagNumber
      report_tmp[1, "Formula"] <- Formula
      report_tmp[1, "Metabolites"] <- Metabolites
      report_tmp[1, "ID"] <- ID
      for(j in 1:sampleNumber){
        alignedGroup_tmp <- alignedGroup[j, ]
        peakGroup <- alignedGroup_tmp$peakGroup[[1]]
        for(m in as.vector(peakGroup$channelIdx)){
          report_tmp[1, paste0(j, "_", m)] <- peakGroup$intb[peakGroup$channelIdx == m]
        }
      }
      report <- rbind(report, report_tmp)
    }
  }
  
  # 5.2 输出report表格
  {
    openxlsx::write.xlsx(report, file = paste0(res_dir, "report.xlsx"))
  }
}
# ****** 5. Result Report ****** #
