program TESSERACT
  !written by Matthew A. Kelly
use number_generators
implicit none
Logical Run, Do_Test_Phase, mewhort_percepts
Integer, Allocatable :: Matrix(:,:),Test(:,:), left(:), indices(:), invIndices(:), &
                        perm(:), permutations(:,:,:), permInput(:), group(:), ranks(:)
Real,  Allocatable :: memory(:,:),percepts(:,:),sentence(:,:),positions(:,:),echo(:),  &
                        orderVec(:),ngrams(:),context(:),placeholder(:),similarities(:), &
                        experience(:),thisProbe(:), probes(:,:), scrambles(:,:),         &
                        units(:,:), unitInput(:)
                        
Character*20, Allocatable ::  Labels(:)
Real time_gauss,time_fft,Time_start,Time_end,Time_sort, SD
Real   result(3)
Character*20 :: name = "minerva"
Character*48 Label_file,    Matrix_file,  Left_file, Group_file,   & ! Input files
             Percepts_file, Place_file, Corpus_file,  & 
             Test_file,     Control_file,  Next_case, &
             scramble_file, memory_file, probe_file, test_matrix, &
             Percepts_out,  Memory_out, Place_out,Concepts_out,   & ! Output files
             Left_out, probe_out, scramble_out, result_out,       &
             matrix_out, test_out, group_out, sent_out, var_out,  &
             unit_file, unit_out, perm_file, perm_out, this_perm_file   ! tesseract files
             
Character*48 :: blank = "                                "
Integer  item,rank,  permsPerSent, numCorrect, sumRank, Vector_Size, Num_Labels
Integer  first,last,cosine_sample
Integer  Inval_count            ! Counts of < 1 word index (short sentences)
Integer  NC                     ! Essentially m choose two
Integer ad1,ad2                 ! used to check array bounds
Integer Da_stat                 ! De-allocate status (currently ignored)
Integer Ios                     ! I/O status - printed and/or ignored
Integer i,j,t,u,h,r,s,p         ! Index
Integer sample                  ! number of 'sample' data to print
Logical :: matrix_only = .False.
Logical :: sparse      = .False.
Logical notLastPerm, even
!                                                
! Control value defaults
!
Integer :: NUM_ITER = 1       ! Number of iterations over which the echo is computed
Integer :: run_num = 0        ! the "run" number, used for output files
Integer :: Num_units = 200    ! total number of units in the tesseract
Integer :: num_unit_files = 1 ! number of files that the tesseract is distributed across
Integer :: POWER = 9          ! MINERVA exponent/power for computing echo
Integer :: NN                 ! dimensionality
Integer :: Max_words = 7      ! sentence length
Integer :: MM                 ! number of unique items in the corpus
Integer :: Num_sent = 125000  ! Sentence count (study set)
Integer :: Num_test = 200     ! Sentence count (test set)
Integer :: Layers =     4     ! Layers default
Integer :: Window_size= 4     ! Window size
Logical :: mewhort_percepts = .false. !if true, transpose percept matrix on read from file
!
Namelist /control/Max_words ,  Vector_Size, cosine_sample,  &
                Percepts_out,  Concepts_out,                &
                Window_size,   Num_sent,    Num_Labels, Layers,     Label_file ,  &
                Percepts_file, Matrix_file, Left_file,  Group_file, Place_file,   &
                Corpus_file,   Next_case,   sample,     matrix_only,sparse,       &
                Memory_out,    Place_out,   test_file,  num_test,                 &
                Left_out,      probe_out,  scramble_out, result_out, matrix_out,  &
                POWER, name, Num_units, unit_file, num_unit_files, unit_out,      &
                perm_file, perm_out, scramble_file, memory_file, probe_file,      &
                test_out, sent_out, test_matrix, run_num, group_out, var_out,     &
                NUM_ITER, mewhort_percepts
!
! Start timing
!
Print*,"MINERVA/tesseract utility"
Print*,"version 'Green Deviations'" 
Print*,"updated Oct. 29th 2018"
Call Cpu_time(Time_start)
Inval_count = 0
Run = .True.
!
! Default file names are blank
!
label_file    = blank; matrix_file  = blank;   left_file    = blank; 
percepts_file = blank; test_file    = blank;   next_case    = blank;
percepts_out  = blank; memory_out   = blank;   place_out    = blank;
place_file    = blank; corpus_file  = blank;   result_out   = blank;
Left_out      = blank; probe_out    = blank;   scramble_out = blank;
group_file    = blank; concepts_out = blank;   matrix_out   = blank;
unit_file     = blank; unit_out     = blank;   perm_file    = blank;
scramble_file = blank; memory_file  = blank;   probe_file   = blank;
perm_out      = blank; test_matrix  = blank;   test_out     = blank;
group_out     = blank; sent_out     = blank;   var_out      = blank;
control_file  = "tesseract.txt"

! initialize random number generator
call RandSeed
SD = Gaussian(0.0, 1.0)					! throw away 1st call to Gaussian

!
!  Clear timing values
!
time_gauss= 0. ; time_fft= 0. ; time_sort = 0.
!
Do while (run)        ! Program is set up to read multiple control files
   Print*, ' Reading control block' 
   Open (7,file=control_file,IOSTAT=Ios)
   Print*,'Ios=',Ios
   Read (7,nml=control,IOSTAT=Ios)
   Close(7) 
   Print*,'Ios=',Ios
   If(Ios < -1 ) then
      Print*,"Bad or missing control block ",control_file
      Print*,"Run terminated"
      Stop
   Endif                                             
   Print*, 'Control block read completed: values merged with defaults: ' 
   Write(6,nml=control)
   If(next_case == blank) then  ! No more control files
      run = .false.             ! Finish this run and then stop
   else
      control_file = next_case
   endif

!  
! Assign control block values
   NN = Vector_size      
   MM = Num_labels
!   
! allocate some stuff
   Allocate (echo(NN*2))
   Allocate (labels(MM))
   Allocate (matrix(Num_sent,Max_words ))
   Allocate (test(Num_test,Max_words))
!
   labels(1) = Blank      !  Used to signal presence/absence of      
                          !  of an input labels file
! 
! READING THE STUDY SET AND TEST SET FROM FILES
!
   Print*, ' Reading label block'
   Call Read_block_CH(labels,label_file,MM,1,sample)
   If(corpus_file .ne. Blank) then   
      Call Read_corpus(matrix,labels,corpus_file,Num_sent,Max_words,MM)
   Else  
      Print*, ' Reading matrix block'
      Call Read_block_IN(matrix,matrix_file,Num_sent,Max_words,sample)
      Print*,'With labels assigned to matrix block sample:'
      Do i = 1, Min(Num_sent,50)   !  Print a sample of the matrix block
        Write(6,'(16A10)') labels(matrix(i,:))
      Enddo
   Endif
   if(matrix_out .ne. Blank) then
        Print*, ' Writing matrix block to file'
        Call Write_block_IN(matrix(1:Num_sent,:),matrix_out,Num_sent,Max_words,run_num)    
   endif
   ! read test set file
    If(test_file .ne. Blank) then      
        Do_Test_Phase = .true.
        ! read corpus needs to be modified to make it work
        Call Read_corpus(test,labels,test_file,Num_test,Max_words,MM)
    elseif(test_matrix .ne. Blank) then
        Do_Test_Phase = .true.
        Print*, ' Reading test matrix file'
        Call Read_block_IN(test,test_matrix,Num_test,Max_words,sample)
        Print*,'With labels assigned to matrix block sample:'
        Do i = 1, Min(Num_test,50)   !  Print a sample of the test matrix
            Write(6,'(16A10)') labels(test(i,:))
        Enddo            
    else
        Do_Test_Phase = .false.
    endif
    If(test_out .ne. Blank) then
        Print*, ' Writing test set matrix to file'
        Call Write_block_IN(test(1:Num_test,:),test_out,Num_test,Max_words,run_num)        
    endif
!
! Allocate a whole bunch of stuff!
!
   Allocate   (left(NN))
   Allocate   (group(NN*2))
   Allocate   (placeholder(NN))
   Allocate   (experience(NN))
   Allocate   (percepts(MM,NN))
   Allocate   (probes(Num_test,NN))
   Allocate   (thisProbe(NN*2))
   Allocate   (orderVec(NN))
   Allocate   (ngrams(NN))   
   Allocate   (memory(Num_sent,NN*2))   
   Allocate   (positions(Max_words,NN))   
   Allocate   (sentence(Max_words,NN))   
   Allocate   (context(NN))
   permsPerSent = factorial(MAX_WORDS)
   Allocate   (scrambles(permsPerSent*Num_test,NN))
   Allocate   (similarities(permsPerSent))
   Allocate   (indices(permsPerSent))
   Allocate   (invIndices(permsPerSent))
   Allocate   (ranks(NUM_TEST))
!
! CONSTRUCT MODEL

!
! create permutations
!  
  If (Left_file == blank) then    
    Print*, 'Random left permutation created' 
    left  = getShuffle(NN)
  Else
    Print*, 'Reading left perm block'          
    Call Read_block_IN(left,left_file,NN,1,sample)
    Print*,'Ios=',Ios
  Endif
  Print*," Maxval(left)=",Maxval(left)
  Print*," Minval(left)=",minval(left)
  Print*," Avg(left)=",Sum(left)/NN
  Print*," Left has ",count((left == 0)),' zero values'

    !  If left_out is non-blank output left
  If(left_out .ne. Blank) then      
        Call Write_block_IN(left,left_out,NN,1,run_num)
  endif

    ! Group permutation (dimensionality = NN * 2)
  If (group_file == blank) then    
    Print*, 'Random group permutation created' 
    group  = getShuffle(NN*2)
  Else
    Print*, 'Reading group perm block'          
    Call Read_block_IN(group,group_file,NN*2,1,sample)
    Print*,'Ios=',Ios
  Endif
  Print*," Maxval(group)=",Maxval(group)
  Print*," Minval(group)=",minval(group)
  Print*," Avg(group)=",Sum(group)/(NN*2)
  Print*," group has ",count((group == 0)),' zero values'

    !  If left_out is non-blank output left
  If(group_out .ne. Blank) then      
        Call Write_block_IN(group,group_out,NN*2,1,run_num)
  endif

! use random vectors as percepts
! these may be generated or read in
!
!
!  Handle initial set of percepts
!  
  If (percepts_file == blank) then
    percepts = 0.
    do i=1,MM
         percepts(i,:) = normalVector(NN) 
    enddo     
    Print*, 'Random percepts block created (Normal Distribution)'
  Else     
    Print*, 'Reading percepts block'
    if (mewhort_percepts) then
        Print*, 'Using Mewhort percepts file ...'
        Call Read_Mewhort(percepts,percepts_file,MM,NN)
        do i=1,MM
            percepts(i,:) = vecNorm(percepts(i,left),NN)
        enddo
        Print*, 'Checking samples from file ...'
        Print*, Labels(120), Labels(523), vectorCosine(percepts(120,:),percepts(523,:),NN) 
        Print*, Labels(120), Labels(536), vectorCosine(percepts(523,:),percepts(536,:),NN) 
        Print*, Labels(523), Labels(536), vectorCosine(percepts(523,:),percepts(536,:),NN)         
    else
        Print*, 'Reading percepts from file ...'
        Call Read_matrix_unf_rl(percepts,percepts_file,MM,NN)
    endif
  Endif
  Print*," Maxval(percepts)=",Maxval(percepts)
  Print*," Minval(percepts)=",minval(percepts)
  Print*," Avg(percepts)=",Sum(percepts)/(MM*NN)
  Print*," Percepts have",count((Percepts == 0)),' zero values' 
  
  !  If Percepts_out is non-blank output percepts
  If(Percepts_out .ne. Blank) then   
    Call Write_unf_rl(percepts,percepts_out,MM*NN,run_num)
  endif
!
!  Handle placeholder
!  
  If (Place_file == blank) then
    placeholder = normalVector(NN)    
    Print*, 'Random placeholders created (Normal Distribution)'
  Else       
    Call Read_block_RL(placeholder,place_file,NN,1,sample)     
    Print*, 'Read completed placeholders block'
  Endif
  Print*," Maxval(placeholder)=",Maxval(placeholder)
  Print*," Minval(placeholder)=",minval(placeholder)
  Print*," Avg(placeholder)=",Sum(placeholder)/NN
  Print*," Placeholder has ",count((placeholder == 0)),' zero values'
  
    !  If Place_out is non-blank output placeholder
  If(place_out .ne. Blank) then      
        Call Write_block(placeholder,place_out,NN,1,run_num)
  endif
  Print*,"placeholder  NN=",NN
  If (NN < 128) stop

! ************************************************************
!           STUDY PHASE
! ************************************************************

    ! build position vectors
    positions(1,:) = placeholder
    do r=2,MAX_WORDS
        positions(r,:) = positions(r-1,left)
    enddo
    Print*,"positions  NN=",NN
    If (NN < 128) stop

    ! read memory from file if file given
    if (memory_file .ne. Blank) then
        Print*, ' Reading memory from file...'         
        Call Read_matrix_unf_RL(memory,memory_file,NUM_SENT,NN*2)
        Print*, ' Finished reading memory from file.'
        Print*," Maxval(memory)=",Maxval(memory)
        Print*," Minval(memory)=",minval(memory)
        Print*," Avg(memory)=",Sum(memory)/(NUM_SENT*NN*2)
        Print*," Memory has",count((memory == 0)),' zero values' 
    else ! otherwise, build memory table
         ! loop over all sentences in the corpus (matrix)
        Print*, ' Building memory from study set...'         
        do s=1,NUM_SENT
            ! build a matrix representing the sentence, one vector per word
            sentence = percepts(matrix(s,:),:) 
            context  = 0. ! build a vector that is a sum of words in sentence
            Do r=1,Max_words
                context = context + sentence(r,:)
            enddo
            ngrams   = hemNgram(sentence,MAX_WORDS,NN,left) ! build a vector of ngrams of words in sentence
            orderVec = getPairs(positions,sentence,MAX_WORDS,NN) ! build a vector of position*sentence pairs
            memory(s,1:NN)       = vecNorm(context,NN)
            memory(s,NN+1:NN+NN) = vecNorm(vecNorm(orderVec,NN) + vecNorm(ngrams,NN),NN)
        enddo        ! end of num_sent loop

        ! write memory table to file
        if (memory_out .ne. Blank) then
            Call Write_unf_rl(memory,memory_out,NUM_SENT*NN*2,run_num)
            Print*,"Memory block written to file ", memory_out
        endif
    endif

! ************************************************************
!   CONSTRUCT TESSERACT
!   (OR READ IT FROM FILE)
! ************************************************************

if (name == 'tessera') then

    Allocate   (units(NUM_UNITS,NN*2))
    Allocate   (permutations(NUM_UNITS,NN*2,POWER))
    Allocate   (permInput(NUM_UNITS*NN*2))
    Allocate   (unitInput(NUM_UNITS*NN*2))

    units        = 0.
    unitInput    = 0.
    permutations = 0
    permInput    = 0

    ! if unit_file and perm_file have been provided, read in tesseract
    if ((unit_file .ne. Blank) .and. (perm_file .ne. Blank)) then
        Print*, ' Reading tesseract from file...'         
        Call Gather_read_rl(unitInput,unit_file,(NUM_UNITS*NN*2)/num_unit_files,num_unit_files)
        units = transpose(reshape(unitInput,(/NN*2,NUM_UNITS/)))
        
        Print*," Maxval(units)=",Maxval(units)
        Print*," Minval(units)=",minval(units)
        Print*," Avg(units)=",Sum(units)/(NUM_UNITS*NN*2)
        Print*," Units have",count((units == 0)),' zero values' 
        
        Print*, ' Reading permutations from file...'         
        do p=1,POWER
            this_perm_file = "P" // Char(Ichar("0")+p) // "_" // perm_file(1:46)
            Call Gather_read_in(permInput,this_perm_file,(NUM_UNITS*NN*2)/num_unit_files,num_unit_files)
            permutations(:,:,p) = transpose(reshape(permInput,(/NN*2,NUM_UNITS/)))
        enddo
        
        Print*," Maxval(permutations)=",Maxval(permutations)
        Print*," Minval(permutations)=",minval(permutations)
        Print*," Avg(permutations)=",Sum(float(permutations))/float(NUM_UNITS*NN*2*POWER)
        Print*," Permutations have",count((permutations == 0)),' zero values' 
        
        Print*, ' Finished reading in units and permutations.'         
    else ! build the tesseract
     
        ! build permutations for the tesseract
        Print*, ' Building permutations for the tesseract...'         
        do u=1,NUM_UNITS
            do p=1,POWER
                permutations(u,:,p) = getShuffle(NN*2)        
            enddo
        enddo
        
        ! add memory traces to the tesseract
        Print*, ' Adding memory traces to the tesseract...' 
        units = 0.
        do s=1,NUM_SENT
            units = tesseractAdd(vecNorm(memory(s,group),NN*2),units,permutations,POWER,NN*2,NUM_UNITS)
        enddo
        Print*, ' Finished adding memory traces to the tesseract' 
    endif
    if (unit_out .ne. Blank) then
        Print*, ' Writing the tesseract to file...' 
        Call Write_unf_rl(transpose(units),unit_out,NUM_UNITS*NN*2,run_num)
        Print*, ' Finished writing the tesseract to file'
    endif
    if (perm_out .ne. Blank) then
        Print*, ' Writing permutations to file...' 
        do p=1,POWER
            this_perm_file = "P" // Char(Ichar("0")+p) // "_" // perm_out(1:46)
            Call Write_unf_IN(transpose(permutations(:,:,p)),this_perm_file,NUM_UNITS*NN*2,run_num)
        enddo
        Print*, ' Finished writing permutations to file' 
    endif    
endif

! ************************************************************
!           TEST PHASE
! ************************************************************

    If(Do_Test_Phase) then
        ! get probes from file if file provided
        if (probe_file .ne. Blank) then
            print*, "Reading probe vectors from file..."
            Call Read_matrix_unf_RL(probes,probe_file,NUM_TEST,NN)
            print*, "Finished reading probes from file."
        else  ! otherwise, generate probes
            print*, "Creating probe vectors..."
            do s=1,NUM_TEST
                 ! build a matrix representing the sentence, one vector per word
                sentence = percepts(test(s,:),:) 
                context  = 0. ! build a vector that is a sum of words in sentence
                Do r=1,Max_words
                    context = context + sentence(r,:)
                enddo
                probes(s,:) = vecNorm(context,NN)
            enddo
            print*, "Finished creating probe vectors."
            ! write probes
            if (probe_out .ne. Blank) then
                Print*,"Writing probe vectors to file..." 
                Call Write_unf_rl(probes,probe_out,NUM_TEST*NN,run_num)
                Print*,"Probe block written" 
            endif
        endif
        
        ! check to see if there's any bad indices in test indices matrix
        do s=1,NUM_TEST
            do r=1,MAX_WORDS
                If (test(s,r) == 0) then
                    Print*,"s,r,test(s,r)=",s,r,test(s,r)
                    NUM_TEST = s - 1
                    Print*,"Zero index detected. Setting NUM_TEST to ",NUM_TEST
                endif
            enddo
        enddo
        
        ! read scrambled sentences from file if a file is given
        if (scramble_file .ne. Blank) then
            print*, "Reading scramble vectors from file..."
            Call Read_matrix_unf_rl(scrambles,scramble_file,permsPerSent*NUM_TEST,NN)
            print*, "Finished reading scramble vectors from file."
        else ! otherwise generate scrambled sentences
            Print*,"Computing scrambled sentence vectors..." 
            t = 1
            do s=1,NUM_TEST
                notLastPerm = .false.
                perm = (/(i,i=1,MAX_WORDS)/)
                do r=1,permsPerSent
                    call nexper(MAX_WORDS,perm,notLastPerm,even)
                    sentence = percepts(test(s,perm),:)
                    ngrams   = hemNgram(sentence,MAX_WORDS,NN,left) ! build a vector of ngrams of words in sentence
                    orderVec = getPairs(positions,sentence,MAX_WORDS,NN) ! build a vector of position*sentence pairs
                    scrambles(t,:) = vecNorm(ngrams,NN) + vecNorm(orderVec,NN)
                    t = t + 1
                enddo
            enddo
        
            ! write scrambled sentences
            If(scramble_out .ne. Blank) then
                Print*,"Writing scrambled sentences to file..."
                Call Write_unf_rl(scrambles,scramble_out,permsPerSent*NUM_TEST*NN,run_num)
                Print*,"Finished writing scrambled sentences."
            Endif            
        endif !end if generate scrambled sentences
        
        print*," " 
        print*,"Testing ",permsPerSent*NUM_TEST," scrambled sentences"
        print*,"generated from ",NUM_TEST," test items"
        print*,"against echoes of ",NN*2," dimensions"
        print*,"generated from ",NUM_SENT," traces in memory"
        print*,"weighted by similarity raised to the power ",POWER
        if (name == 'tessera') then
            print*,"stored across ",NUM_UNITS," distributed units"
            print*,"utilizing ",NUM_UNITS*POWER," permutations"
        endif
        print*," "
        
        ! generate echoes and test against scrambled sentences
        print*, "Computing echoes..."
        t = 1
        numCorrect = 0
        sumRank    = 0
        
        ! open file for writing sentences generated to file
        if (sent_out .ne. blank) then
            Open(9,file=sent_out,IOSTAT=Ios)
        endif
        if (var_out .ne. blank) then
            Open(8,file=var_out,IOSTAT=Ios)
        endif
        ! begin testing
        do s=1,NUM_TEST
            ! pad probe with zeroes
            thisProbe = 0.
            thisProbe(1:NN) = probes(s,:)

            ! create echo by probing with context (unordered words)
            ! create echo
            do i=1,NUM_ITER
                if (name == 'minerva') then
                    if (POWER == 0) then
                        echo = beagleEcho(thisProbe,memory,NN*2,NUM_SENT)
                    else
                        echo = minervaEcho(thisProbe,memory,POWER,NN*2,NUM_SENT)
                    endif
                elseif (name == 'tessera') then
                    echo(group) = tesseractEcho(vecNorm(thisProbe(group),NN*2),units,permutations,POWER,NN*2,NUM_UNITS)
                else
                    Print*, 'No model specified: Error in echo generation.'
                    Stop
                endif
                thisProbe(NN+1:NN+NN) = vecNorm(echo(NN+1:NN+NN),NN)
            enddo

            ! test latter half of echo (order) against each scrambled sentence
            do r=1,permsPerSent
                similarities(r) = vectorCosine(echo(NN+1:NN*2),scrambles(t,:),NN)
                t = t + 1
            enddo
            ! print cosine of the correct answer
            Print*,"     "
            Print*,"sim (1, unsorted)=",similarities(1)
            if (var_out .ne. blank) then
                ! write normalized deviation to file
                write(8,'(F30.20)') deviation(similarities(1),similarities)
            endif
            ! sort cosines
            Call Sort_index(similarities,indices,permsPerSent,permsPerSent,'dn')
            ! print cosine of the highest ranked answer
            Print*,"sim (",indices(1),", sorted)=",similarities(indices(1))            
            invIndices(indices) = (/(i,i=1,permsPerSent)/)
            rank = invIndices(1)
            Print*,"rank ",rank
            if (rank==1) then
                numCorrect = numCorrect + 1
            endif

            notLastPerm = .false.
            perm = (/(i,i=1,MAX_WORDS)/)
            do i=1,indices(1)
                call nexper(MAX_WORDS,perm,notLastPerm,even)
            enddo
            print*,labels(test(s,perm))
            
            if (sent_out .ne. blank) then     
                write(9,'(I4, 2x, 7A20)') rank, labels(test(s,perm))
            endif
            ranks(s) = rank
            sumRank  = sumRank + rank
        enddo ! loop over test items 1 to NUM_TEST
        if (sent_out .ne. blank) then
            close(9) ! close sentence file
        endif
        if (var_out .ne. blank) then
            close(8) ! close variance file
        endif
        
        Print*,"Finished testing echoes and getting results" 
        result(1) = float(numCorrect) / float(NUM_TEST) ! percentage correct
        result(2) = float(sumRank)    / float(NUM_TEST) ! average ranking
        result(3) = median(float(ranks),NUM_TEST)       ! median ranking
        
        !write result
        Print*,"percentage correct = ",result(1)
        Print*,"mean rank = ",result(2)
        Print*,"median rank = ",result(3)
        If(result_out .ne. Blank) then
            Call Write_block(result,result_out,3,1,run_num)
        endif
    endif ! endif Do_Test_Phase
enddo ! do while run

! END MAIN PROGRAM
! BEGIN FUNCTION STATEMENTS
CONTAINS
!
!***************************************************************
!
SUBROUTINE NEXPER(N, A, MTC, EVEN)
INTEGER N
INTEGER, DIMENSION(N) :: A
INTEGER S, D, NM3, IA, I1, L, I, J, M
LOGICAL MTC, EVEN
    IF (MTC) GOTO 10
    NM3 = N-3
    DO 1 I=1,N
1   A(I)=I
    MTC=.TRUE.
5   EVEN=.TRUE.
    IF(N .EQ. 1) GOTO 8
6   IF(A(N) .NE. 1 .OR. A(1) .NE. 2+MOD(N,2)) RETURN
    IF(N .LE. 3) GOTO 8
    DO 7 I=1,NM3
    IF(A(I+1) .NE. A(I)+1) RETURN
7   CONTINUE
8   MTC=.FALSE.
    RETURN
10  IF(N .EQ. 1) GOTO 27
    IF(.NOT. EVEN) GOTO 20
    IA=A(1)
    A(1)=A(2)
    A(2)=IA
    EVEN=.FALSE.
    GOTO 6
20  S=0
    DO 26 I1=2,N
25  IA=A(I1)
    I=I1-1
    D=0
    DO 30 J=1,I
30  IF(A(J) .GT. IA) D=D+1
    S=D+S
    IF(D .NE. I*MOD(S,2)) GOTO 35
26  CONTINUE
27  A(1)=0
    GOTO 8
35  M=MOD(S+1,2)*(N+1)
    DO 40 J=1,I
    IF(ISIGN(1,A(J)-IA) .EQ. ISIGN(1,A(J)-M)) GOTO 40
    M=A(J)
    L=J
40  CONTINUE
    A(L)=IA
    A(I1)=M
    EVEN=.TRUE.
    RETURN
END SUBROUTINE NEXPER
!
!***************************************************************
!
real function median(dataset,size)
    implicit none
    integer, intent(in)                  :: size
    real,    intent(in), dimension(size) :: dataset
    integer,             dimension(size) :: indices
    logical :: odd
    integer :: middle
    
    Call Sort_index(dataset,indices,size,size,'up')
    odd    = (mod(size,2) == 1)
    middle = size / 2
    if (odd) then
        median = dataset(indices(middle+1))
    else
        median = dataset(indices(middle))
        median = median + dataset(indices(middle+1))
        median = median / 2.
    endif
end function median
!
!***************************************************************
!
integer function factorial(n)
implicit none
integer, intent(in) :: n
integer :: i, ans

ans = 1
do i = 1,n
    ans = ans * i
enddo
factorial = ans

end function factorial
!
!***************************************************************
!
Function cconv(V1,V2,Nfeatures)
    implicit none
    integer, intent(in)                         :: Nfeatures
    real, dimension(Nfeatures), intent(in)      :: V1, V2
    real, dimension(Nfeatures)                  :: cconv

    cconv = ifft(fft(V1,Nfeatures) * fft(V2,Nfeatures),Nfeatures)
End function cconv
!
!***********************************************************************
!
Function fft(a,n)
    implicit none
    include 'fftw3.f'
    integer,            intent(in) :: n
    real, dimension(n), intent(in) :: a
    complex, dimension(n)          :: a_in, a_out, fft
    integer                        :: plan
    !print*,'in fft'
    a_in = dcmplx(real(a), 0.)
    
    !print*,'executing fft'
    call dfftw_plan_dft_1d(plan,n,a_in,a_out,FFTW_FORWARD,FFTW_ESTIMATE)
    call dfftw_execute_dft(plan,  a_in,a_out)
    call dfftw_destroy_plan(plan)
    !print*,'finished fft'
    fft = a_out
    !print*,'leaving fft'
End function fft
!
!***********************************************************************
!
Function ifft(a_in,n)
    implicit none
    include 'fftw3.f'
    integer,               intent(in)        :: n
    complex, dimension(n), intent(in)        :: a_in
    complex, dimension(n)                    :: a_out 
    real,           dimension(n)             :: ifft
    integer                                  :: plan
    
    call dfftw_plan_dft_1d(plan,n,a_in,a_out,FFTW_BACKWARD,FFTW_ESTIMATE)
    call dfftw_execute_dft(plan,  a_in, a_out)
    call dfftw_destroy_plan(plan)
    
    ifft = real(a_out, 8)
    ifft = ifft/n
End function ifft
!
!*************************************************************************
!
function getShuffle( n )
use number_generators
!
!  GETSHUFFLE provides a permutation and inverse permutation 
!             of the numbers 1:N.
!   perm: for shuffling a thing, shuffledThing = thing(perm)
!   invPerm: for unshuffling, thing = shuffledThing(invPerm)
!
   Integer :: getShuffle (n)
   Real    :: vector(n)
   Integer n,i
   
   ! generate n real values from 0...1
   do i=1,n
        vector(i) = flat(1.) 
   enddo
   
   Call Sort_index(vector,getShuffle,n,n,'dn') !(~,perm) = sort(rand((1,N)))        
                                               !(~,invPerm) = sort(perm);
end function getShuffle 
!
!***************************************************************
!
Function hemNgram( window, WL, NN, leftin )
!hemNgram produces all the convolutional n-grams of row vectors in percepts
!   n-grams vary in size from 1 to numOfItems   
!   Contrast with hemUOG, which gets all skip-grams / unconstrained open
!   grams of row vectors in percepts.
!
!   window: matrix of environmental vectors of dimensions [numOfItems,N]
!             where: numOfItems is the number of vectors
!                    N is the dimensionality of each vector
!
!   chunk: sum of all n-grams of the environmental vectors in window
!          
!   left: permutation for indicating that a vector is the left operand of
!         a non-commutative circular convolution. By default, the function
!         uses no permutation (i.e., the numbers 1 to N in ascending order)
!
!   WL: window length
!   NN: vector dimensionality

Integer NN,WW,WL,i
Real     :: hemNgram(NN)
Real     :: gram(NN)
Real     :: window(WL,NN) 
Integer             :: left(NN)
Integer,optional    :: leftin(NN)

! check if a permutation has been passed in
if( present(leftin) )then
    left = leftin
else ! if not, do not permute (i.e., use 1 to the dimensionality of the vectors)
    left = (/(i,i=1,NN)/)
endif

gram      = window(1,:)
hemNgram  = window(1,:)

Do i=2,WL
    gram     = window(i,:) + cconv(gram(left),window(i,:),NN)
    hemNgram = hemNgram + gram
enddo

end function hemNgram
!
!***************************************************************
!
Function getPairs(setA, setB, WL, NN)
! computes all pairs of setA(i,:) * setB(i,:) and sums
Integer NN,WW,WL,i,p
Real     :: getPairs(NN)
Real     :: setA(WL,NN) 
Real     :: setB(WL,NN) 

getPairs = 0.
Do i=1,WL
    getPairs = getPairs + cconv(setA(i,:),setB(i,:),NN)
enddo

end function getPairs
!
!***************************************************************
!
Integer Function Lookup(Table,Value,Nvalues,mode)
    Character*20 Table(Nvalues),Value, revalue,table_entry
    Character*6  mode
    Integer Nvalues, i, L1,L2
    revalue = adjustl(Value)
    L1 = index(revalue,' ')
    revalue(L1:20) = Blank
    Lookup = 0
    Do i = 1,Nvalues
       table_entry = adjustl(Table(i))
       L1 = index(table_entry,' ')
       If(L1 > 0) table_entry(L1:20) = Blank
       If(table_entry == reValue) then
          Lookup = i
          Return
       Endif
    Enddo
    If(mode=="create") then
       Nvalues = Nvalues + 1
       Table(Nvalues) = Value
       Lookup = Nvalues
       Return
    Endif 
End function lookup
!
!***************************************************************
!
Function Lbl(indx)
    Character*10 Lbl
    Integer      indx
    Lbl = Blank
    If(indx > 0 ) Lbl = labels(indx)
End function Lbl

!
!***************************************************************
!   MINERVA ECHO retrieves an ECHO from the MINERVA MEMORY MODEL
!   an ECHO is an aggregate representation of stored memories
!   a weighted sum of all vectors/traces in memory
!   each weighted by their similarity to the vector PROBE
!   each similarity raised to an exponent POWER
!   TRACE_DIM is the dimensionality of each memory trace
!   NUM_TRACES is the number of memory traces in memory
!***************************************************************
!      
function minervaEcho(probe,memory,POWER,TRACE_DIM,NUM_TRACES)
integer i, TRACE_DIM, NUM_TRACES, POWER
real    activation, probe(TRACE_DIM), minervaEcho(TRACE_DIM), &
	    memory(NUM_TRACES,TRACE_DIM) 
	    
minervaEcho = 0.
do i=1,NUM_TRACES
    activation  = vectorCosine(probe,memory(i,:),TRACE_DIM)**POWER
    minervaEcho = minervaEcho + activation * memory(i,:)
enddo

end function minervaEcho
!
!***************************************************************
!   BEAGLE ECHO finds the memory trace vector most similar  
!   to the PROBE in MEMORY
!   Each memory trace is TRACE_DIM dimensions
!   There are NUM_TRACES memory traces    
!***************************************************************
!
function beagleEcho(probe,memory,TRACE_DIM,NUM_TRACES)
integer i, TRACE_DIM, NUM_TRACES, best
real    similarity, highestSim, probe(TRACE_DIM), beagleEcho(TRACE_DIM), &
	    memory(NUM_TRACES,TRACE_DIM) 

highestSim = 0.
do i=1,NUM_TRACES
    similarity  = vectorCosine(probe,memory(i,:),TRACE_DIM)
    if (similarity > highestSim) then
        highestSim = similarity
        best = i
    endif
enddo

beagleEcho = memory(best,:)

end function beagleEcho
!
!***************************************************************
!
function tesseractAdd(x,units,permutations,POWER,TRACE_DIM,NUM_UNITS)
    Implicit none
    integer, intent(in) :: TRACE_DIM, NUM_UNITS, POWER, &
                           permutations(NUM_UNITS,TRACE_DIM,POWER)
    real, intent(in)    :: x(TRACE_DIM), units(NUM_UNITS,TRACE_DIM)
    real                :: tesseractAdd(NUM_UNITS,TRACE_DIM)
    complex, dimension(TRACE_DIM) :: trace
    integer u,p

    ! add trace to each units unit
    do u=1,NUM_UNITS    
        ! compute circular convolution:
        ! trace = x * P1 x * P2 x * P3 x
        trace = fft(x,TRACE_DIM)
        do p=1,POWER
            trace = trace * fft(x(permutations(u,:,p)),TRACE_DIM)
        enddo
        tesseractAdd(u,:) = units(u,:) + ifft(trace,TRACE_DIM)
    enddo
end function tesseractAdd
!
!***************************************************************
!
function tesseractEcho(probe,units,permutations,POWER,TRACE_DIM,NUM_UNITS)
    Implicit none
    integer, intent(in) :: TRACE_DIM, NUM_UNITS, POWER, &
                           permutations(NUM_UNITS,TRACE_DIM,POWER)
    real, intent(in)    :: probe(TRACE_DIM), units(NUM_UNITS,TRACE_DIM)
    integer             :: i, u, p, inv(TRACE_DIM)

    real,             dimension(TRACE_DIM) :: echo, tesseractEcho
    double precision, dimension(TRACE_DIM) :: echoSum
    complex,   dimension(TRACE_DIM)        :: powerProbe
    
    ! inv is the permutation (1, N ... 2) where N is the dimensionality
    inv = (/1,(i,i=TRACE_DIM,2,-1)/)
    ! average over all units of the model
    !print*,'initialization of echo sum'
    echoSum = 0.
    do u=1,NUM_UNITS
        ! compute F P1 Pinv probe
        !print*,'initialization of power probe'
        powerProbe = 1.
        !print*,'building power probe'
        do p=1,POWER
            ! compute F P2...j Pinv probe
            powerProbe = powerProbe * fft(probe(permutations(u,inv,p)),TRACE_DIM)
        enddo
        !print*,'computing echo'
        ! compute the echo
        echo = ifft(powerProbe * fft(units(u,:),TRACE_DIM),TRACE_DIM)
        !print*,'summing across echoes'
        ! sum across all echoes
        echoSum = echoSum + echo
    enddo
    !print*,'taking the mean of the echoes'
    ! take the mean of the echoes
    tesseractEcho = echoSum / NUM_UNITS
end function tesseractEcho
!
!***************************************************************
!
Function normalVector( n )
use number_generators
!
! Generate a vector of n random values
! selected from an normal distribution
! with a mean of 0 and a variance of 1/n
! (i.e. a standard deviation of sqrt(1/n)).
! The vector is then normalized to have a mean of zero
! and Euclidean length (i.e magnitude) of one.
! mean and standard deviation of the normal distribution
!
Real   :: normalVector(n)
Real      vector(n)
Real   mn,sd,De_norm
Real :: start, finish
Integer i,n
Call cpu_time(start)  
sd = sqrt(1./n)
mn = 0.
!
! sample n random values from the normal distribution to create a vector
!
Do i=1,n
    Vector(i) = gaussian(mn,sd)
enddo
!
! The vector's mean and magnitude should be respectively close to
! 0 and 1 given the properties of the distribution from which the
! values of the vector have been sampled. Nonetheless, here we
! ensure that the mean is exactly 0 and magnitude is exactly 1.
!
mn = sum(vector)/n
vector = vector - mn  
De_norm = sqrt(dot_product(vector,vector))
If ( De_norm > 2. .or. De_norm < .5) then
     Print*,'De_norm =',De_norm
     Print '(20f8.4)',vector
     stop
endif
normalVector = vector / sqrt(dot_product(vector,vector))
Call Cpu_time(finish)
Time_gauss = Time_gauss + finish - start
end function normalVector
!
!
!**************************************************************************
!
Subroutine Read_corpus(Mtrx,Label,File_name,NS,NW,MM)
!
! read a text corpus - return a matrix (Mtrx(i,j)
! where the value of Mtrx(i,j) is the index in Label
! of the j'th word in the i'th sentence of the corpus 
Integer NS                     !  Number of sentences
Integer NW                     !  Word count this sentence
Integer MM                     !  Number of labels ( words)
Integer Mtrx(NS,NW)            !  Output matrix
Integer         i              !  Sentence index
Integer         max_i          !  Index of last sentence found
Integer         j              !  Word index current line
Integer         rli            !  Word index in create mode
Integer      oldrli            !  Prior valu of rli
Integer         lbi            !  Word index in label
Integer   ::    sti = 0        !  sentence index (count)
Integer         nerr           !  error count in lookup
Character*48 File_name
Character*20  Label(MM)
Character*20  Sent (NW)        !  Sentence as read
Integer      iSent (NW)        !  Sentence as an index
Integer ::    Stats(NW)
Integer max_j,max_max_j
Character*6 Mode
!
Mode = "exists"                         ! Mode for lookup is governed by
If(Label(1)== Blank)Mode = "create"     ! whether we are creating the label
                                        ! file or using an existing one
nerr   = 0; max_i = 0;   max_j =0; 
oldrli = 0; max_max_j=0; sti   =0
rli = MM
If(Mode == "create") then
  rli = 0
  Label = Blank
Endif
Stats =0
Open (8,file=File_name,IOSTAT=Ios)
Print*,'Read_corpus on formatted open Ios=',Ios,' on file ',File_name
Print*,"Read_corpus attempting to acquire=",NS," sentences."
Print*,"Each sentence has maximum  length:",NW
Print*,"Unique words list (labels) length:",MM
!
Do i = 1,NS   !   We are trying to get NS sentences from the corpus
   Sent = Blank
   Read (8,*,IOSTAT=Ios)Sent
   if (Ios < 0) then
        exit
   else
       Do j = 1,NW   !   j'th word in the sentence
          If (Sent(j).ne. Blank) then
              lbi  = Lookup(Label,Sent(j),rli,mode)
              if (lbi == 0) then        
                    nerr = nerr + 1
                    If(nerr.lt.20) print*,"error lookup: not found ",Sent(j)
                    If(nerr > 40) then
                        Print*,"Program terminated: max not founds is 40"  
                        Stop
                    endif
              endif
              iSent(j) = lbi
              max_i = i
              max_j = j
          Else
              iSent(j) = 0
          Endif
       Enddo
       If(.not.Sparse) then
           if (i == max_i) mtrx(max_i,:) = iSent
           max_max_j = max(max_max_j,max_j)
           stats(max_j) = stats(max_j) + 1
       else                                  ! Sparse is true and the
         If(rli > oldrli ) then              ! sentence has new words
            sti = sti + 1                    ! increment sentence count
            mtrx(sti,:) = iSent              ! otherwise not saved
            max_max_j = max(max_max_j,max_j) ! when Sparse is set
            stats(max_j) = stats(max_j) + 1
            oldrli = rli
         endif ! rli > oldrli
       endif ! sparse / not sparse
    endif ! ios < 0
Enddo
NS = max_i
MM = rli
Print*,"Read_corpus  sentences read     = ",NS 
If(Sparse) then
	Print*,"Sparse is set. Sentences cut    = ",NS - sti
	NS = sti
	Print*,"Sentence count has been reset   = ",NS
Endif 
Print*,"Unique words list re-set to length",MM
Print*,"Sentence statistics: length, number of occurrences"
Print'(8(I4,I6,4x))',(j,Stats(j),j=1,max_max_j)
Print*,'Sample of file:'
!Do i = 1, Min(NS,50)   !
Do i = Max(1,NS - 50), NS
   Write(6,'(12A10)') (Lbl(Mtrx(i,j)),j=1,Min(NW,12))
Enddo
Close(8)
Print*,"Read_corpus completed"
End Subroutine Read_corpus
!
!
!***************************************************************************
!
Subroutine  Sort_index(x,index,lindx,sample,direction)
!
! This subroutine receives an array x() and creates an index into 
! it such that x(index(i)) is in in direction order, i=1,sample.
      Implicit  none
      Integer                       :: lindx    ! array length
      Integer                       :: sample   ! length of req'd sample
      Integer                       :: i,j,loc  ! array indices
      Integer                       :: NumPairs ! progress counter
      Integer                       :: LastSwap ! index of last swap
      Integer                       :: Temp     ! index temporary store
      Real  ,   dimension(lindx)    :: x        ! data being indexed
      Real  ,   dimension(lindx)    :: xx       ! a copy of x  
      Integer,  dimension(lindx)    :: index    ! index to data times
      Character*2 direction         !  sort may be up or down
      Real start, finish
!
      Call Cpu_time(Start)
!     
      index =(/(i,i=1,lindx)/)
!
! Decreasing size  e.g. biggest  first is selected by direction = "dn"
! Increasing size  e.g. smallest first is selected by direction = "up"
! To handle up we just flips the signs before and after
!
      If(direction == "up") x = -x
      If (Sample < Lindx/20 ) then
         xx = x
         Do i = 1, sample 			  
            Loc = maxloc(xx,1)   ! find max 
            index(i)= Loc        ! save the index of the maximum
            xx(Loc) = -Huge(0.0) ! destroy the prior maximum in xx
         End do 
      Else
      NumPairs = lindx - 1
      Index =(/(i,i=1,lindx)/)  ! Initialze index to 1,2,3,4…..

     Do
     If (NumPairs == 0) Exit

     ! If no more pairs to check, terminate repetition

     ! Otherwise scan the sublist of the first NumPairs pairs in 
     ! the list, interchanging items that are out of order

     LastSwap = 1
     Do I = 1, NumPairs
        If (x(Index(i)) < x(Index(i+1))) THEN
           ! Items out of order -- interchange them
           Temp = Index(i)
           Index(i) = Index(i+1)
           Index(i+1) = Temp         
           LastSwap = i      ! Record position of last swap
        Endif
     Enddo
     NumPairs = LastSwap - 1
  Enddo
  Endif
  If(direction == "up") x = -x 
  Call Cpu_time(Finish)
  Time_sort = Time_sort + Finish -Start
!
End subroutine  Sort_index
!
!**************************************************************************
!
function vecNorm(a, NN)
!
! vecNorm normalizes vector "a" to a Euclidean length (i.e., magnitude) of 1
!
    Integer NN
    Real    :: vecNorm(NN)
    Real    a(NN), aa
!
    aa = dot_product(a,a)
    if (aa == 0) aa = 1
    vecNorm  = a / sqrt(aa)
end function vecNorm
!
!**************************************************************************
!
Real   function  vectorCosine(x, y, Np)
    Integer Np
    Real   x(Np),y(Np),lengthX,lengthY,productXY
!
! VECTOR COSINE computes the cosine of the angle between the vectors x & y
! and outputs a value between 1 and -1.
! The cosine is useful as a measure of similarity:
! 0 means the vector are orthogonal, or completely dissimilar
! +1 means the vectors are identical
! -1 means the vectors are exact opposites
! Special case: If one of the vectors has a magnitude of zero, the cosine is zero.
    lengthX = dot_product(x,x) 
    lengthY = dot_product(y,y) 
    productXY = dot_product(x,y) 
    if(lengthX * lengthY == 0 )then 
        vectorCosine = 0;
    else
        vectorCosine = productXY / (sqrt(lengthX) * sqrt(lengthY))
    endif
end function vectorCosine

! Computes the sample standard deviation
real function std(vector)
      implicit none
      ! input
      real, dimension(:), intent(in) :: vector
      ! internal
      integer :: n
      real    :: mean, variance

      n         = size(vector)
      mean      = sum(vector) / float(n)
      variance  = sum((vector - mean)**2) / float(n-1)
      std       = sqrt(variance)
end function std

! Computes the deviation of a data point
real function deviation(datum, vector)
    implicit none
    ! input
    real,               intent(in) :: datum
    real, dimension(:), intent(in) :: vector
    ! internal
    integer :: n
    real    :: mean

    n         = size(vector)
    mean      = sum(vector) / float(n)
    deviation = (datum - mean) / std(vector)
end function deviation


!
!
End Program
!******************************************************************
!******************************************************************
!******************************************************************
!
!  The following code is deliberately not in the CONTAINS block
!
!
Subroutine Read_block_IN(Array,File_name,L1,L2,sample)
!
! read a block of INTEGER values
!
Integer i,j,L1,L2,sample
Integer Array(L1,L2)
Character*48 File_name
If(Index(File_name,"unf") > 0 .or. Index(File_name,"UNF") > 0 ) then 
        Open (8,file=File_name,Form='unformatted',IOSTAT=Ios)
        Print*,'Unformatted Open Ios=',Ios,' on file ',File_name
        Read (8,iostat=Ios)Array
else 
        Open (8,file=File_name,IOSTAT=Ios)
        Print*,'Formatted open Ios=',Ios,' on file ',File_name
        Read (8,*,IOSTAT=Ios)Array
        Array= reshape(Array,(/L1,L2/),order =(/2,1/))
Endif
Print*,'On read of ',file_name,' Ios=',Ios
Print*," Maxval ",File_name,Maxval(Array)
Print*," Minval ",File_name,Minval(Array)
If(sample > 0) then
 Print*,'Sample of file:'
 If(L2 > 1) then
   Do i = 1, Min(L1,12)    
      Write(6,'(12I8)') (Array(i,j),j=1,Min(L2,sample))
   Enddo
 Else    
      Write(6,'(12I8)') (Array(i,1),i=1,Min(L1,sample))
 Endif
Endif
Close(8)
End Subroutine Read_block_IN
!
Subroutine Read_block_RL(Array,File_name,L1,L2,sample)
!
! read a block of Real   values
!
Integer i,j,L1,L2,sample
Real   Array(L1,L2)
Character*48 File_name
IF(File_name == "                                ") Return
If(Index(File_name,"unf") > 0 .or. Index(File_name,"UNF") > 0 ) then 
        Open (8,file=File_name,Form='unformatted',IOSTAT=Ios)
        Print*,'Unformatted Open Ios=',Ios,' on file ',File_name
        Read (8,iostat=Ios)Array
else 
        Open (8,file=File_name,IOSTAT=Ios)
        Print*,'Formatted Open Ios=',Ios,' on file ',File_name
        Read (8,*,IOSTAT=Ios)Array
!        Array= reshape(Array,(/L1,L2/),order =(/2,1/))
Endif
Print*,'On read of ',file_name,' Ios=',Ios
If (sample > 0 ) then
 Print*,'Sample of file:'
 If(L2 > 1) then
    Do i = 1, Min(L1,12)   !
       Write(6,'(2x,12F10.5)') (Array(i,j),j=1,Min(L2,sample))
    Enddo
 Else
    Write(6,'(2x,12F10.5)') (Array(i,1),i=1,Min(L1,sample))
 Endif
Endif  
Close(8)
End Subroutine Read_block_RL
!
!
Subroutine Read_block_CH(Array,File_name,L1,L2,sample)
!
! read a block of CHAR*20 values
!
Integer i,j,L1,L2,sample
Character*20 Array(L1,L2)
Character*48 File_name
IF(File_name == "                                ") Return
If(Index(File_name,"unf") > 0 .or. Index(File_name,"UNF") > 0 ) then 
        Open (8,file=File_name,Form='unformatted',IOSTAT=Ios)
        Print*,'Unformatted Open Ios=',Ios,' on file ',File_name
        Read (8,iostat=Ios)Array
else 
        Open (8,file=File_name,IOSTAT=Ios)
        Print*,'Formatted Open Ios=',Ios,' on file ',File_name
        Read (8,*,IOSTAT=Ios)Array
        !Array= reshape(Array,(/L1,L2/),order =(/2,1/))
Endif
Print*,'On read of ',file_name,' Ios=',Ios
Close(8)
If (sample > 0 ) then
 Print*,'Sample of file:'
 If ( L2 > 1) then
      Do i = 1, Min(L1,12)   
         Write(6,'(12(1x,A10))') (Array(i,j),j=1,Min(L2,sample))
      Enddo
 Else   
         Write(6,'(12(1x,A10))') (Array(i,1),i=1,Min(L1,sample))
 Endif
Endif
End Subroutine Read_block_CH
!
!
Subroutine Write_block(Array,File_name,L1,L2,layer)
    Real,         intent(in) :: Array(L1,L2)
    Character*48, intent(in) :: File_name
    Integer,      intent(in) :: L1,L2,layer
    Real Copy(L2,L1)
    Character*51 File_name2

    If(File_name == "                                ") Return
    If(layer < 10) then
        File_name2 = "L" // Char(Ichar("0")+layer) // "_" // File_name
    else
        File_name2 = "L" // Char(Ichar("A") + layer - 10) // "_" // File_name
    endif
    If(Index(File_name,"unf") > 0 .or. Index(File_name,"UNF") > 0 ) then 
            Open (11,file=File_name,Form='unformatted',IOSTAT=Ios)
            Print*,'Unformatted Open Ios=',Ios,' on file ',File_name
            Write (11,iostat=Ios)Array
    else 
            Open (11,file=File_name2,IOSTAT=Ios)
            Print*,'Formatted Open Ios=',Ios,' on file ',File_name
            Copy = transpose(Array)
            Write (11,'(20F11.4)',IOSTAT=Ios)Copy
    Endif
    Print*,'On write to ',File_name2,' Ios=',Ios
    Close(11)
End Subroutine Write_block
!
!
Subroutine Write_block_IN(Array,File_name,L1,L2,layer)
    Integer, intent(in)      :: L1,L2,layer,Array(L1,L2)
    Character*48, intent(in) :: File_name    
    Integer                  :: Copy(L2,L1)
    Character*48 format_string
    Character*51 File_name2
    
    If(File_name == "                                ") Return
    File_name2 = "L" // Char(Ichar("0")+layer) // "_" // File_name
    Open (11,file=File_name2,IOSTAT=Ios)
    Print*,'Formatted Open Ios=',Ios,' on file ',File_name
    Copy = transpose(Array)
    Write (format_string,'("(",I2,"I7)")') L2
    Write (11,format_string,IOSTAT=Ios) Copy
    Print*,'On write to ',File_name2,' Ios=',Ios
    Close(11)
End Subroutine Write_block_IN

! GATHER READ SUBROUTINES
! _rl is for real*4
! _r8 is for real*8
! _in is for integer, calls rl

Subroutine Gather_read_rl(Array,File_name,L1,Nfiles)
Real   Array(L1*Nfiles)
Real   Recrd(L1)
Integer:: L1,Nfiles
Character*48 File_name
Character*55 File_name2
If(File_name == "                                ")then
   Print*,"No file name specified."
   Return
Endif
J = 1
Do I = 1,Nfiles
   Write(File_name2,"(I2,A,A)")I,File_name(1:len_trim(File_name)),".unf"
   File_name2 = adjustl(File_name2) 
   Open (11,file=File_name2,Form='unformatted',IOSTAT=Ios)
   !Print*,'Unformatted Open Ios=',Ios,' on file ',file_name2,Ios
   If(Ios < 0) then
     Print*,"Error on open of file:",File_name2,Ios
     Close(11)
     Return
   Endif
   Read (11,iostat=Ios)Recrd
   If(Ios < 0) then
      Print*,"Error on read of file:",File_name2,Ios
      Close(11)
      return
   Endif
   Array(j:(j + L1 - 1)) = Recrd
   J = J + L1
   Close(11)
Enddo
End Subroutine Gather_read_rl
!
!
Subroutine Gather_read_in(Array,File_name,L1,Nfiles)
    implicit none
    Integer    Array(L1*Nfiles)
    Integer :: L1,Nfiles,i,j,ios
    Character*48 File_name
    Integer   Recrd(L1)
    Character*55 File_name2

    If(File_name == "                                ")then
       Print*,"No file name specified."
       Return
    Endif
    J = 1
    Do I = 1,Nfiles
       Write(File_name2,"(I2,A,A)")I,File_name(1:len_trim(File_name)),".unf"
       File_name2 = adjustl(File_name2) 
       Open (11,file=File_name2,Form='unformatted',IOSTAT=Ios)
       !Print*,'Unformatted Open Ios=',Ios,' on file ',file_name2,Ios
       If(Ios .ne. 0) then
         Print*,"Error on open of file:",File_name2,Ios
         Close(11)
         Return
       Endif
       Read (11,iostat=Ios)Recrd
       If(Ios .ne. 0) then
          Print*,"Error on read of file:",File_name2,Ios
          Close(11)
          return
       Endif
       Array(j:(j + L1 - 1)) = Recrd
       J = J + L1
       Close(11)
    Enddo
End Subroutine Gather_read_in


! Unformatted write subroutines
!   _rl is for real*4
!   _r8 is for real*8
!   _in is for integer (calls _rl)

Subroutine Write_unf_rl(Array,File_name,L1,Indx)
    Real   Array(L1)
    Integer:: L1,Indx
    Character*48 File_name
    Character*55 File_name2
    If(File_name == "                                ")then
       Print*,"No file name specified."
       Return
    Endif
    Write(File_name2,"(I2,A,A)")Indx,File_name(1:len_trim(File_name)),".unf"
    File_name2 = adjustl(File_name2) 
    Print*,File_name2
    Open (11,file=File_name2,Form='unformatted',IOSTAT=Ios)
    !Print*,'Unformatted Open Ios=',Ios,' on file ',file_name2,Ios
    Write (11,iostat=Ios)Array
    If(Ios < 0) then
       Print*,"Error on write to file:",File_name2,Ios
       Stop
    Endif
    Close(11)
End Subroutine Write_unf_rl
!
!
Subroutine Write_unf_IN(Array,File_name,L1,Indx)
    Integer Array(L1)
    Integer:: L1,Indx
    Character*48 File_name
    Character*55 File_name2
    If(File_name == "                                ")then
       Print*,"No file name specified."
       Return
    Endif
    Write(File_name2,"(I2,A,A)")Indx,File_name(1:len_trim(File_name)),".unf"
    File_name2 = adjustl(File_name2) 
    Print*,File_name2
    Open (11,file=File_name2,Form='unformatted',IOSTAT=Ios)
    !Print*,'Unformatted Open Ios=',Ios,' on file ',file_name2,Ios
    Write (11,iostat=Ios)Array
    If(Ios < 0) then
       Print*,"Error on write to file:",File_name2,Ios
       Stop
    Endif
    Close(11)
End Subroutine Write_unf_IN
!
!
Subroutine Read_matrix_unf_rl(Array,File_name,L1,L2)
    Implicit none
    Real,         intent(out) :: Array(L1,L2)
    Character*48, intent(in)  :: File_name
    Integer,      intent(in)  :: L1,L2
    Integer Ios

    If(File_name == "                                ")then
       Print*,"No file name specified."
       Return
    Endif
    Open (11,file=File_name,Form='unformatted',IOSTAT=Ios)
    If(Ios .ne. 0) then
        Print*,"Error on open of file:",File_name,Ios
        Close(11)
        Return
    Endif
    Read (11,iostat=Ios) Array
    If(Ios .ne. 0) then
          Print*,"Error on read of file:",File_name,Ios
          Close(11)
          return
    Endif
    Close(11)
End Subroutine Read_matrix_unf_rl

Subroutine Read_matrix_unf_IN(Array,File_name,L1,L2)
    Implicit none
    Integer,      intent(out) :: Array(L1,L2)
    Character*48, intent(in)  :: File_name
    Integer,      intent(in)  :: L1,L2
    Integer Ios

    If(File_name == "                                ")then
       Print*,"No file name specified."
       Return
    Endif
    Open (11,file=File_name,Form='unformatted',IOSTAT=Ios)
    If(Ios .ne. 0) then
        Print*,"Error on open of file:",File_name,Ios
        Close(11)
        Return
    Endif
    Read (11,iostat=Ios) Array
    If(Ios .ne. 0) then
          Print*,"Error on read of file:",File_name,Ios
          Close(11)
          return
    Endif
    Close(11)
End Subroutine Read_matrix_unf_IN

! Read Mewhort's BEAGLE visual / environmental vectors in as percepts
Subroutine Read_Mewhort(percepts,File_name,MM,NN)
    Implicit none
    Real,         intent(out) :: percepts(MM,NN)
    Character*48, intent(in)  :: File_name
    Integer,      intent(in)  :: MM,NN
    Integer i
    
	percepts = 0.0
	open(1, FILE = File_name, status = 'old', form = 'unformatted')
    write(*, *) ' Read percept vectors'  
    do i = 1, MM
            read(1) percepts(i, :)
    enddo
    close(1)
End Subroutine