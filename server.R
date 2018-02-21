library(shiny)
library(shinythemes)
library(DT)
library(ggplot2)
library(gridExtra)



source('PowerCalc_Rare.r')
source('calc_whole.r')
source('PowerCalc_RareSample_A.R')

cutrange <- function(start,end,size){
  len <- end-start
  dif <- len/(size-1)
  result <- rep(0,size)
  for(i in 1:size){
    result[i] <- start+(i-1)*dif
  }
  return(result)
}

grid_transform = function(grid){
  if(grid=="Fast"){
    return(20)
  }else if(grid=="Intermediate"){
    return(50)
  }else{
    return(100)
  }
}


# Define UI for slider demo application

function(input, output) {
  
  # Reactive expression to compose a data frame containing all of
  # the values
  
  # Compose data frame
  
  
  #,input$nameEsseble,input$JJ
  
  
  
  data1 <- eventReactive(input$update,{
    if(input$TypeofCalculation=='Sample Size Calculation'){
      if(input$SNPoption_s=="Single SNP"){
        if(input$evoption_s=='Single EV'){
          result <- get_Aprox_Sample(EV=(input$EV_s)/100,
                           PowerThr=input$PowerThreshold_s,
                           alpha = input$Alpha_s,
                           ONESNP=T)

          
          
          
          
        }else{
          EV_new <- cutrange(input$EV_range_s[1],input$EV_range_s[2],5)
          
          get_Aprox_Sample(EV=(EV_new)/100,
                           PowerThr=input$PowerThreshold_s,
                           alpha = input$Alpha_s,
                           ONESNP=T)
        }
        
      }else{
        if(input$evoption_s=='Single EV'){
          if((input$method_s)!='Burden Test'){
            
            
            get_Aprox_Sample(EV=(input$EV_s)/100,
                             PowerThr=input$PowerThreshold_s,
                             alpha = input$Alpha_s,
                             TEST=input$method_s,
                             QT='CC',
                             PC=input$PC_s,
                             JJ=input$JJ_s,
                             PRC=input$PRC_s)
            
            
            
          }
          else if((input$method_s)=='Burden Test'){
            
            
            get_Aprox_Sample(EV=(input$EV_s)/100,
                             PowerThr=input$PowerThreshold_s,
                             alpha = input$Alpha_s,
                             TEST=input$method_s,
                             PC=input$PC_new_s,
                             JJ=input$JJ_s,
                             PRC=input$PRC_s)
            
            
          }
          
        }else if(input$evoption_s=='Range of EV'){
          if(input$method_s!='Burden Test'){
            EV_new <- cutrange(input$EV_range_s[1],input$EV_range_s[2],5)
         
            
            get_Aprox_Sample(EV=(EV_new)/100,
                             PowerThr=input$PowerThreshold_s,
                             alpha = input$Alpha_s,
                             TEST=input$method_s,
                             PC=input$PC_s,
                             JJ=input$JJ_s)
            
          }
          else if(input$method_s=='Burden Test'){
            EV_new <- cutrange(input$EV_range_s[1],input$EV_range_s[2],5)
            get_Aprox_Sample(EV=(EV_new)/100,
                             PowerThr=input$PowerThreshold_s,
                             alpha = input$Alpha_s,
                             TEST=input$method_s,
                             PC=input$PC_new_s,
                             JJ=input$JJ_s,
                             PRC=input$PRC_s)
            
          }
        }
        
        
      }
      
    }else{
      if(input$SNPoption=="Single SNP"){
        if(input$evoption=='Single EV'){
          get_Aprox((input$EV)/100,
                    input$Alpha,
                    (input$ncases+input$ncont),
                    input$ncases,input$ncont,
                    ONESNP=T)
        }else{
          EV_new <- cutrange(input$EV_range[1],input$EV_range[2],5)
        
          get_Aprox((EV_new)/100,
                    input$Alpha,
                    (input$ncases+input$ncont),
                    input$ncases,input$ncont,
                    ONESNP=T)
        }
        
      }else if(input$SNPoption=='Whole Genome'){
        
        calcGenomeLevel(K=input$K,
                        m=input$m,
                        grid=grid_transform(input$grid),
                        epr=(input$GEV/(input$K*100)),
                        alphaT=input$Alpha_whole,
                        Total=as.numeric((input$ncases+input$ncont)),
                        CASE=as.numeric(input$ncases),
                        CONTROL=as.numeric(input$ncont),
                        TEST = input$method_whole,
                        QT='CC',
                        PC=input$PC_whole,
                        JJ=input$JJ_whole)
        
      }else{
        if(input$evoption=='Single EV'){
          if((input$method)!='Burden Test'){
            
            get_Aprox(EV=(input$EV)/100,
                      alpha=input$Alpha,
                      Total=(input$ncases+input$ncont),
                      CASE=input$ncases,
                      CONTROL=input$ncont,
                      TEST=input$method,
                      QT='CC',
                      PC=input$PC,
                      JJ=input$JJ,
                      PRC=input$PRC)
            
            
            
            
            
          }
          else if((input$method)=='Burden Test'){
            get_Aprox(EV=(input$EV)/100,
                      alpha=input$Alpha,
                      Total=(input$ncases+input$ncont),
                      CASE=input$ncases,
                      CONTROL=input$ncont,
                      PC=input$PC_new,
                      TEST=input$method, JJ=input$JJ,PRC=input$PRC)
            
          }
          
        }else if(input$evoption=='Range of EV'){
          if(input$method!='Burden Test'){
            EV_new <- cutrange(input$EV_range[1],input$EV_range[2],5)
            
            get_Aprox(EV=(EV_new)/100,
                      alpha=input$Alpha,
                      Total=(input$ncases+input$ncont),
                      CASE=input$ncases,
                      CONTROL=input$ncont,
                      TEST=input$method,
                      PC=input$PC,
                      JJ=input$JJ)
          }
          else if(input$method=='Burden Test'){
            EV_new <- cutrange(input$EV_range[1],input$EV_range[2],5)
            
            get_Aprox(EV=(EV_new)/100,
                      alpha=input$Alpha,
                      Total=(input$ncases+input$ncont),
                      CASE=input$ncases,
                      CONTROL=input$ncont,
                      TEST=input$method,
                      PC=input$PC_new,
                      JJ=input$JJ,
                      PRC=input$PRC)
          }
        }
        
        
      }
      
    }
    
  })
  
  
  
  
  
  
  data2 <- eventReactive(input$update2,{
    if(input$TypeofCalculation2=='Sample Size Calculation'){
      if(input$SNPoption2_s=="Single SNP"){
        if(input$evoption2_s=='Single EV'){
          get_Aprox_Sample(EV=(input$EV2_s)/100,
                           PowerThr=input$PowerThreshold2_s,
                           alpha = input$Alpha2_s,
                           QT='QT',
                           ONESNP=T)
        }else{
          EV_new2 <- cutrange(input$EV_range2_s[1],input$EV_range2_s[2],5)
          get_Aprox_Sample(EV=(EV_new2)/100,
                           PowerThr=input$PowerThreshold2_s,
                           alpha = input$Alpha2_s,
                           ONESNP=T,
                           QT='QT')
        }
        
      }else{
        if(input$evoption2_s=='Single EV'){
          if((input$method2_s)!='Burden Test'){
            get_Aprox_Sample(EV=(input$EV2_s)/100,
                             PowerThr=input$PowerThreshold2_s,
                             alpha = input$Alpha2_s,
                             TEST=input$method2_s,
                             QT='QT',
                             PC=input$PC2_s,
                             JJ=input$JJ2_s,
                             PRC=input$PRC2_s)
            
            
          }
          else if((input$method2_s)=='Burden Test'){
            get_Aprox_Sample(EV=(input$EV2_s)/100,
                             PowerThr=input$PowerThreshold2_s,
                             alpha = input$Alpha2_s,
                             TEST=input$method2_s,
                             PC=input$PC_new2_s,
                             JJ=input$JJ2_s,
                             PRC=input$PRC2_s,
                             QT = 'QT')
            
            
          }
          
        }else if(input$evoption2_s=='Range of EV'){
          if(input$method2_s!='Burden Test'){
            EV_new2 <- cutrange(input$EV_range2_s[1],input$EV_range2_s[2],5)
           
            get_Aprox_Sample(EV=(EV_new2)/100,
                             PowerThr=input$PowerThreshold2_s,
                             alpha = input$Alpha2_s,
                             TEST=input$method2_s,
                             PC=input$PC2_s,
                             JJ=input$JJ2_s,
                             QT='QT')
            
          }
          else if(input$method2_s=='Burden Test'){
            EV_new2 <- cutrange(input$EV_range2_s[1],input$EV_range2_s[2],5)
            
            get_Aprox_Sample(EV=(EV_new2)/100,
                             PowerThr=input$PowerThreshold2_s,
                             alpha = input$Alpha2_s,
                             TEST=input$method2_s,
                             PC=input$PC_new2_s,
                             JJ=input$JJ2_s,
                             PRC=input$PRC2_s,
                             QT='QT')
          
          }
        }
        
        
      }
      
      
      
      
    }else{
      if(input$SNPoption2=="Single SNP"){
        if(input$evoption2=='Single EV'){
          get_Aprox((input$EV2)/100,
                    input$Alpha2,
                    (input$total2),
                    QT='QT',
                    ONESNP=T)
        }else{
          EV_new2 <- cutrange(input$EV_range2[1],input$EV_range2[2],5)
          
          get_Aprox((EV_new2)/100,
                    input$Alpha2,
                    (input$total2),
                    QT='QT',
                    ONESNP=T)
        }
        
      }else if(input$SNPoption2=='Whole Genome'){
        calcGenomeLevel(K=input$K2,
                        m=input$m2,
                        grid=grid_transform(input$grid2),
                        epr=(input$GEV2/(input$K2*100)),
                        alphaT=input$Alpha_whole2,
                        Total=as.numeric(input$total2),
                        TEST = input$method_whole2,
                        QT='QT',
                        PC=input$PC_whole2,
                        JJ=input$JJ_whole2)
      }else{
        if(input$evoption2=='Single EV'){
          if((input$method2)!='Burden Test'){
            
            get_Aprox(EV=(input$EV2)/100,
                      alpha=input$Alpha2,
                      Total=input$total2,
                      PC=input$PC2,
                      TEST=input$method2, 
                      JJ=input$JJ2,
                      PRC=input$PRC2,
                      QT='QT')
            
          }
          else if((input$method2)=='Burden Test'){
            get_Aprox(EV=(input$EV2)/100,
                      alpha=input$Alpha2,
                      Total=(input$total2),
                      PC=input$PC_new2,
                      TEST=input$method2, 
                      JJ=input$JJ2,
                      PRC=input$PRC2,
                      QT='QT')
            
          }
          
        }else if(input$evoption2=='Range of EV'){
          if(input$method2!='Burden Test'){
            EV_new2 <- cutrange(input$EV_range2[1],input$EV_range2[2],5)
            
            
            get_Aprox(EV=(EV_new2)/100,
                      alpha=input$Alpha2,
                      Total=input$total2,
                      TEST=input$method2,
                      PC=input$PC2,
                      JJ=input$JJ2,
                      QT='QT')
          }
          else if(input$method2=='Burden Test'){
            EV_new2 <- cutrange(input$EV_range2[1],input$EV_range2[2],5)
            
            
            get_Aprox(EV=(EV_new2)/100,
                      alpha=input$Alpha2,
                      Total=input$total2,
                      TEST=input$method2,
                      PC=input$PC_new2,
                      JJ=input$JJ2,
                      PRC=input$PRC2,
                      QT='QT')
          }
        }
        
        
      }
      
      
      
      
    }
    
    
    
  })
  
  
  
  
  
  
  
  busyIndicator <- function (text = "Loading..", img = "busyIndicator/ajaxloaderq.gif", wait = 1000) {
    tagList(singleton(tags$head(tags$link(rel = "stylesheet", 
                                          type = "text/css", href = "busyIndicator/busyIndicator.css"))), 
            div(class = "shinysky-busy-indicator", p(text), img(src = img)), 
            tags$script(sprintf("\tsetInterval(function(){\n  \t\t \t if ($('html').hasClass('shiny-busy')) {\n  \t\t      
                                setTimeout(function() {\n  \t\t      if ($('html').hasClass('shiny-busy')) {\n  \t\t        $
                                ('div.shinysky-busy-indicator').show()\n  \t\t      }\n  \t\t    }, %d)  \t\t    \n  \t\t  } 
                                else {\n  \t\t    $('div.shinysky-busy-indicator').hide()\n  \t\t  }\n  \t\t},100)\n  \t\t", 
                                wait)))
    }
  
  
  output$values <- DT::renderDataTable( 
    
    withProgress( message = "WE ARE COMPUTING..", value = 0.9,{ 
      data1()[[1]]
    })  
    
    ,options=list(dom='t',autoWidth=F,columnDefs = list(list(width = "150", targets = "_all"))),rownames=F
  )
  output$plot = renderPlot(
    grid.arrange(data1()[[2]],data1()[[3]],
                 data1()[[4]],ncol=1)
    ,
    height=700,
    width=566
  )
  
  
  
  
  
  output$values2 <- DT::renderDataTable( 
    
    withProgress( message = "WE ARE COMPUTING..", value = 0.9,{ 
      data2()[[1]]
    })  
    
    ,options=list(dom='t',autoWidth=F,columnDefs = list(list(width = "150", targets = "_all"))),rownames=F
  )
  output$plot2 = renderPlot(
    grid.arrange(data2()[[2]],data2()[[3]],
                 data2()[[4]],ncol=1)
,
height=700,
width=566
    )
  
  # output$values2_r <- DT::renderDataTable( 
  #   withProgress(message="We are computing based on different EV range",{ 
  #     data2_r()[[1]]
  #   })  
  #   ,options=list(dom='t',autoWidth=T,columnDefs = list(list(width = '150px', targets = "_all"))),rownames=F
  # )
  # output$plot2_r <- renderPlot(data2_r()[[2]],height = 300,width = 600)
  #  
  
  
  
  
  }