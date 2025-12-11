# BIBLIOTECAS
library(R.matlab)

####################################################################################################
# FUNCOES 

# as funcoes mais importantes sao mostradas abaixo, funcoes menores serao escritas no proprio script

  # Separando condições com duas partes (rest_1 e 2 para HBN e task_pre e pos para INPD)
separating_conditions_interval <- function(file){
  
  condition1 <- list()
  condition2 <- list()
  
  # rest_1 e rest_2 HBN
  if (is.list(file)) { 
    
    for (i in 1:length(file[[1]])) {
      half <- length(file[[1]][[i]]) / 2
      
      # lidando com matrizes multicolunas (white matter e CSF)
      if (length(file[[1]][[i]]) > length(file[[1]][[1]])) {
        
        condition1[[i]] <- as.matrix(file[[1]][[i]][1:(length(file[[1]][[1]])/2), ], 1)
        condition2[[i]] <- as.matrix(file[[1]][[i]][((length(file[[1]][[1]])/2) + 1):length(file[[1]][[1]]), ], 1)
      } 
      
      # restante dos dados
      else {
        
        condition1[[i]] <- as.matrix(file[[1]][[i]][1:half, ], 1)
        condition2[[i]] <- as.matrix(file[[1]][[i]][(half + 1):length(file[[1]][[i]]), ], 1)
      }
    }
  }
  
  else {
    
    # task_pre e task_pos INPD
    matfile <- readMat(file)
    
    for (j in 1:length(matfile$data)){
      half_matfile <- length(matfile$data[[j]][[1]]) / 2
      
      # lidando com matrizes multicolunas (white matter e CSF)
      if (length(matfile$data[[j]][[1]]) > length(matfile$data[[1]][[1]])) {
        
        condition1[[j]] <- as.matrix(matfile$data[[j]][[1]][1:(length(matfile$data[[1]][[1]])/2), ], 1)
        condition2[[j]] <- as.matrix(matfile$data[[j]][[1]][((length(matfile$data[[1]][[1]])/2) + 1):length(matfile$data[[1]][[1]]), ], 1)
      } 
      
      # restante dos dados
      else {
        
        condition1[[j]] <- as.matrix(matfile$data[[j]][[1]][1:half_matfile, ], 1)
        condition2[[j]] <- as.matrix(matfile$data[[j]][[1]][(half_matfile + 1):length(matfile$data[[j]][[1]]), ], 1)
      }
      
      
    }
  }
  
  condition1 <- list(list(condition1))
  condition2 <- list(list(condition2))
  
  return(c(condition1, condition2))
}

  # separando condições especificas ('movieDM', 'movieTP' e 'peer', todas do HBN) 
reduzir_data <- function(matlab_file, condition) {
  
  matfile <- readMat(matlab_file)
  data <- matfile$data
  conditionweights <- matfile$conditionweights
  
  # define índices válidos (>0)
  cw <- conditionweights[[condition]][[1]]
  idx <- which(cw > 0)
  
  # função auxiliar para aplicar idx corretamente
  data_reduced <- vector("list", length(data))
  for (i in seq_along(data)) {
    x <- data[[i]]
    if (is.list(x) && length(x) == 1) x <- x[[1]]
    if (!is.null(x) && is.matrix(x)) {
      n <- nrow(x)
      valid_idx <- idx[idx <= n]
      if (length(valid_idx) > 0) {
        # mantém todas as colunas, selecionando apenas as linhas válidas
        data_reduced[[i]] <- x[valid_idx, , drop = FALSE]
      }
    }
  }
  
  return(list(data_reduced))
}

  # metaestabilidade adaptada para condições em partes e condicoes especificas
metastability_adapted <- function(matriz, inicio_atlas=168, final_atlas=500, tamanho_janela=NULL) {
  
  
  parcellation_matrix <- as.data.frame(matrix(nrow=nrow(as.data.frame(matriz[[1]][[1]])), ncol=length(matriz[[1]])))
  
  for (i in 1:length(matriz[[1]])) {
    x <- matriz[[1]][[i]]
    if (is.matrix(x) && ncol(x) == 1) { # matrizes com somente uma coluna
      parcellation_matrix[, i] <- x
    } else {
      parcellation_matrix[, i] <- NA  # casos com mais de uma coluna
    }
  }
  
  gordon_matrix <- parcellation_matrix[, inicio_atlas:final_atlas]
  
  gordon_matrix <- gordon_matrix[, !names(gordon_matrix) %in% names(gordon_matrix)[colSums(gordon_matrix == 0, na.rm=TRUE) == nrow(gordon_matrix) & colSums(is.na(gordon_matrix)) == 0]]
  
  # cálculo da coerência espacial
  mean_fmri <- c()
  diff_fmri <- matrix(NA, nrow=nrow(gordon_matrix), ncol=ncol(gordon_matrix))
  V <- c()
  
  # média p/ cada TR (linha)
  for (i in 1:nrow(gordon_matrix)) {
    mean_fmri[i] <- mean(as.numeric(gordon_matrix[i, ]))
  }
  
  # matriz dos desvios e coerência espacial
  for (i in 1:nrow(gordon_matrix)) {
    for (j in 1:ncol(gordon_matrix)) {
      diff_fmri[i, j] <- abs(gordon_matrix[i, j] - mean_fmri[i])
    }
    V[i] <- sum(diff_fmri[i, ]) # somando os desvios de cada coluna vs media para cada linha
  }
  V <- V/ncol(gordon_matrix) 
  
  
  # ajustando tamanho da janela caso não seja passado como argumento
  if (is.null(tamanho_janela)) {
    tamanho_janela <- nrow(gordon_matrix)
  }
  
  # metaestabilidade
  V <- as.data.frame(V) 
  number_windows <- floor(nrow(gordon_matrix)/tamanho_janela)
  metastability <- matrix(NA, nrow=number_windows, 1)
  
  start <- 1
  end <- tamanho_janela
  for (r in 1:number_windows) {
    metastability[r, 1] <- sd(V[start:end, 1], na.rm=TRUE)
    start <- start + tamanho_janela
    end <- end + tamanho_janela
  }
  
  metastability <- as.data.frame(metastability)
  colnames(metastability) <- 'metastability'
  
  return(metastability) 
}

  # aplicando as funcoes 'reduzir_data' e 'metastability_adapted' nos arquivos dos sujeitos
lista_resultados_adaptado <- function(lista, cond) {
  extraindo_condicoes <- list()
  resultado <- list()
  
  for (i in 1:length(lista)) {
    extraindo_condicoes[[i]] <- reduzir_data(lista[[i]], cond)
  }
  
  for (j in 1:length(extraindo_condicoes)) {
    resultado[[j]] <- metastability_adapted(extraindo_condicoes[[j]])
  }
  return(data.frame(resultado))
}


####################################################################################################
# IMPORTACAO

  # INPD
INPD_subjects <- list.files('F:/TASK_INDP_patients/task_INDP/results/preprocessing/', pattern='ROI.*\\Condition000.mat$', full.names=TRUE, recursive=FALSE)

  # HBN

    # rest (condicao em duas partes, separadas na secao de resultados)

      # rest CBIC
CBIC_rest_subjects <- list.files('F:/HBN/CBIC/CBIC/results/preprocessing/', pattern='ROI.*\\Condition001', full.names=TRUE, recursive=FALSE)
CBIC_rest_subjects <- CBIC_rest_subjects[2:3] # sujeito 1 nao tem YBOCS

      # rest RU
RU_rest_subjects <- list.files('F:/HBN/RU/RU/results/preprocessing/', pattern='ROI.*\\Condition001', full.names=TRUE, recursive=FALSE)

    # condicoes especificas (movieDM, movieTP e peer)

      # movieDM  (condicao 2 no arquivo output do CONN)
CBIC_movieDM <- list.files('F:/HBN/CBIC/CBIC/results/preprocessing/', pattern='ROI.*\\Condition002.mat$', full.names=TRUE, recursive=FALSE)
RU_movieDM <- list.files('F:/HBN/RU/RU/results/preprocessing/', pattern='ROI.*\\Condition002.mat$', full.names=TRUE, recursive=FALSE)
movieDM <- c(CBIC_movieDM, RU_movieDM)
movieDM <- movieDM[2:8] # sub_1 sem YBOCS

      # movieTP (condicao 3 no arquivo output do CONN)
CBIC_movieTP <- list.files('F:/HBN/CBIC/CBIC/results/preprocessing/', pattern='ROI.*\\Condition003.mat$', full.names=TRUE, recursive=FALSE)
RU_movieTP <- list.files('F:/HBN/RU/RU/results/preprocessing/', pattern='ROI.*\\Condition003.mat$', full.names=TRUE, recursive=FALSE)
movieTP <- c(CBIC_movieTP, RU_movieTP)
movieTP <- movieTP[2:8] # sub_1 sem YBOCS
movieTP[[8]] <- NA # O elemento 8 (sub5 - condition003) não possui dados funcionais de movieTP

      # peer (condicao 4 no arquivo output do CONN)
CBIC_peer <- list.files('F:/HBN/CBIC/CBIC/results/preprocessing/', pattern='ROI.*\\Condition004.mat$', full.names=TRUE, recursive=FALSE)
RU_peer <- list.files('F:/HBN/RU/RU/results/preprocessing/', pattern='ROI.*\\Condition004.mat$', full.names=TRUE, recursive=FALSE)
peer <- c(CBIC_peer, RU_peer)
peer <- peer[2:8] # sub_1 sem YBOCS

####################################################################################################
# RESULTADOS
  
# metaestabilidade para condicoes em partes (task pre e pos e rest_1 e 2)

  # task_pre e task_pos INPD
resultados_INPD <- c()

for (i in 1:length(INPD_subjects)){
  resultados_INPD[[i]] <- separating_conditions_interval(INPD_subjects[[i]])
}

nomes <- list.dirs('F:/TASK_INDP_patients/task_patients_data/', recursive=FALSE, full.names=FALSE)
names(resultados_INPD) <- nomes

    # metaestabilidade task 1 (pre) INPD
metaestabilidade_task_pre_INPD <- list()
for (i in 1:length(resultados_INPD)){
  metaestabilidade_task_pre_INPD[[i]] <- metastability_adapted(resultados_INPD[[i]][[1]])
}

metaestabilidade_task_pre_INPD <- t(data.frame(metaestabilidade_task_pre_INPD))
colnames(metaestabilidade_task_pre_INPD) <- 'metastability'
metaestabilidade_task_pre_INPD <- data.frame('Sujeitos'=paste0('task_pre_', nomes), metaestabilidade_task_pre_INPD)


    # metaestabilidade task 2 (pos) INPD
metaestabilidade_task_pos_INPD <- list()
for (i in 1:length(resultados_INPD)){
  metaestabilidade_task_pos_INPD[[i]] <- metastability_adapted(resultados_INPD[[i]][[2]])
}

metaestabilidade_task_pos_INPD<- t(data.frame(metaestabilidade_task_pos_INPD))
colnames(metaestabilidade_task_pos_INPD) <- 'metastability'
metaestabilidade_task_pos_INPD <- data.frame('Sujeitos'=paste0('task_pos_', nomes), metaestabilidade_task_pos_INPD)


  # rest_1 e rest_2 CBIC
extracting_rest_CBIC <- list()

    # extraindo o a condição rest do arquivo output do CONN
for (i in 1:length(CBIC_rest_subjects)){
  extracting_rest_CBIC[[i]] <- reduzir_data(CBIC_rest_subjects[[i]], 1) 
}


    # Separando rest1 e rest2 para cada sujeito
results_separating_interval_rest_CBIC <- list()

for (j in 1:length(results_reduzing_data_CBIC_rest)){
  results_separating_interval_rest_CBIC[[j]] <- separating_conditions_interval(extracting_rest_CBIC[[j]])
}


    # nomeando os sujeitos na lista final
names(results_separating_interval_rest_CBIC) <- c('sub_1', 'sub_2')
names(results_separating_interval_rest_CBIC[[1]]) <- c('rest_1', 'rest_2')
names(results_separating_interval_rest_CBIC[[2]]) <- c('rest_1', 'rest_2')


    # calculando metaestabilidade com a lista de rests separados

      # rest_1 CBIC
metaestabilidade_rest1_CBIC <- list()
for (i in 1:length(results_separating_interval_rest_CBIC)){
  metaestabilidade_rest1_CBIC[[i]] <- metastability_adapted(results_separating_interval_rest_CBIC[[i]][[1]])
}

metaestabilidade_rest1_CBIC <- t(data.frame(metaestabilidade_rest1_CBIC))
colnames(metaestabilidade_rest1_CBIC) <- 'Metastability'
metaestabilidade_rest1_CBIC <- data.frame('Sujeito'=c('sub_1', 'sub_2'), metaestabilidade_rest1_CBIC)

      # rest_2 CBIC
metaestabilidade_rest2_CBIC <- list()
for (i in 1:length(results_separating_interval_rest_CBIC)){
  metaestabilidade_rest2_CBIC[[i]] <- metastability_adapted(results_separating_interval_rest_CBIC[[i]][[2]])
}

metaestabilidade_rest2_CBIC <- t(data.frame(metaestabilidade_rest2_CBIC))
colnames(metaestabilidade_rest2_CBIC) <- 'Metastability'
metaestabilidade_rest2_CBIC <- data.frame('Sujeito'=c('sub_1', 'sub_2'), metaestabilidade_rest2_CBIC)


  # rest_1 e rest_2 RU
extracting_rest_RU <- list()

    # extraindo a condicao rest do arquivo output do CONN
for (i in 1:length(RU_rest_subjects)) {
  extracting_rest_RU[[i]] <- reduzir_data(RU_rest_subjects[[i]], 1) 
}

    # separando rest_1 e rest_2
results_separating_interval_rest_RU <- list()

for (j in 1:length(extracting_rest_RU)){
  results_separating_interval_rest_RU[[j]] <- separating_conditions_interval(extracting_rest_RU[[j]])
}


    # nomeando as listas
names(results_separating_interval_rest_RU) <- c('sub_3', 'sub_4', 'sub_5', 'sub_6', 'sub_7')

    # calculado metaestabilidade 

      # rest_1 RU
metaestabilidade_rest1_RU <- list()

for (i in 1:length(results_separating_interval_rest_RU)){
  metaestabilidade_rest1_RU[[i]] <- metastability_adapted(results_separating_interval_rest_RU[[i]][[1]]) 
}

metaestabilidade_rest1_RU <- t(data.frame(metaestabilidade_rest1_RU))
colnames(metaestabilidade_rest1_RU) <- 'Metastability'
metaestabilidade_rest1_RU <- data.frame('Sujeito'=paste0('rest_1_', c('sub_3', 'sub_4', 'sub_5', 'sub_6', 'sub_7')), metaestabilidade_rest1_RU)

      # rest_2 RU
metaestabilidade_rest2_RU <- list()

for (i in 1:length(results_separating_interval_rest_RU)){
  metaestabilidade_rest2_RU[[i]] <- metastability_adapted(results_separating_interval_rest_RU[[i]][[2]]) 
}

metaestabilidade_rest2_RU <- t(data.frame(metaestabilidade_rest2_RU))
colnames(metaestabilidade_rest2_RU) <- 'Metastability'
metaestabilidade_rest2_RU <- data.frame('Sujeito'=paste0('rest_2_', c('sub_3', 'sub_4', 'sub_5', 'sub_6', 'sub_7')), metaestabilidade_rest2_RU)


# metaestabilidde para condicoes especificas (movieDM, movieTP e peer)

results_movieDM <- lista_resultados_adaptado(movieDM, 2) 

results_movieTP <- lista_resultados_adaptado(movieTP, 3)
results_movieTP[[8]] <- NA

results_peer <- lista_resultados_adaptado(peer, 4)

resultados_HBN_especificos <- t(data.frame(c(results_movieDM, results_movieTP, results_peer)))
colnames(resultados_HBN_especificos) <- 'Metastability'
resultados_HBN_especificos <- data.frame(Sujeito=c('movieDM_sub_1', 'movieDM_sub_2','movieDM_sub_3', 'movieDM_sub_4','movieDM_sub_5', 
                                                         'movieDM_sub_6', 'movieDM_sub_7', 'movieTP_sub_1', 'movieTP_sub_2', 'movieTP_sub_3', 
                                                         'movieTP_sub_4', 'movieTP_sub_5', 'movieTP_sub_6', 'movieTP_sub_7', 'peer_sub_1', 
                                                         'peer_sub_2', 'peer_sub_3', 'peer_sub_4', 'peer_sub_5', 'peer_sub_6', 'peer_sub_7'), 
                                         resultados_HBN_especificos)


# juntando dataframes
INPD_results <- rbind(metaestabilidade_task_pre_INPD, metaestabilidade_task_pos_INPD)
HBN_results_em_partes <- rbind(metaestabilidade_rest1_CBIC,
                               metaestabilidade_rest1_RU,
                               metaestabilidade_rest2_CBIC,
                               metaestabilidade_rest2_RU)

names(HBN_results_em_partes) <- names(resultados_HBN_especificos) # evitar erros de nomeacao no rbind
HBN_results <- rbind(HBN_results_em_partes, resultados_HBN_especificos)

names(INPD_results) <- names(HBN_results) 

# dataframe final dos resultados preliminares
resultados_preliminares <- rbind(INPD_results, HBN_results)
