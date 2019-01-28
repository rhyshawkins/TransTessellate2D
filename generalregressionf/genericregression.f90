module genericregression

  !
  ! Calling C/C++ code from fortran and visa-versa can have compiler dependent issues.
  ! This should work when using GNU g++ and gfortran for relatively recent compilers
  ! (tested using 8.2.1).
  !
  use iso_c_binding
  implicit none

  !
  ! The interface expects the user to store their own data. Here we create a structure
  ! to store each observation consisting of an x and y coordinate and value with error
  ! sigma.
  !
  type :: observation_t

     real(kind = c_double), dimension(1) :: x
     real(kind = c_double), dimension(1) :: y
     real(kind = c_double) :: value
     real(kind = c_double) :: sigma

  end type observation_t

  !
  ! The actual observations are stored in a module array
  !
  integer :: nobservations
  type(observation_t), dimension(:), allocatable :: observations


contains

  !
  ! This function is needed for interoperability between C++ and Fortran
  ! to convert and array of character to a string, see the open call in
  ! gvs2_loaddata_
  !
  function copystring(a) 
    character, intent(in) :: a(:)
    character(len = size(a)) :: copystring
    integer :: i
    do i = 1, size(a)
       copystring(i:i) = a(i)
    end do
  end function copystring

  !
  ! This function is called on initialisation and is used to inform the
  ! framework how many 2D models are required and how many hierarchical
  ! parameters. 
  !
  function gvcart_initialise_(nmodels, nhierarchical) bind(c)

    integer, intent(out) :: nmodels
    integer, intent(out) :: nhierarchical

    integer :: gvcart_initialise_

    nmodels = 1
    nhierarchical = 1

    gvcart_initialise_ = 0
    
  end function gvcart_initialise_

  !
  ! The data file to load is passed via command line to this function and here
  ! the user supplied function must load the actual data and inform the
  ! framework of the points required for each observation. This is
  ! done by registering the points with the callback function addobs
  !
  function gvcart_loaddata_(n, filename, addobs) bind(c)

    integer, intent(in) :: n
    character, dimension(n), intent(in) :: filename
    
     interface iaddobs
       function addobs(npoints, modelindex, lons, lats) bind(c)
         use iso_c_binding
         integer :: addobs
         integer(kind = c_int), intent(in) :: npoints
         integer(kind = c_int), intent(in) :: modelindex
         real(kind = c_double), dimension(npoints), intent(in) :: lons
         real(kind = c_double), dimension(npoints), intent(in) :: lats
       end function addobs
    end interface iaddobs
   
    integer :: gvcart_loaddata_

    !
    ! Local parameters
    !
    integer, parameter :: FUNIT = 10
    
    integer :: status
    integer :: i
    
    open(FUNIT, file = copystring(filename), status = 'old', iostat = status)
    if (status .ne. 0) then
       gvcart_loaddata_ = -1
       return
    end if
    
    !
    ! Custom text file format, first line contains the number of points
    ! Each subsequent line has floating point lon, lat, value, sigma
    !
    read(FUNIT, *) nobservations
    
    !
    ! Allocate storage
    !
    allocate(observations(nobservations))

    !
    ! Read each observation in
    !
    do i = 1, nobservations
       
       read(FUNIT, *) observations(i)%x(1), observations(i)%y(1), observations(i)%value, observations(i)%sigma

       !
       ! Here we call addobs to register this observation with the inversion
       ! C++ code. For regression there is one point per observation and
       ! these are stored in the observation_t structure with length 1 arrays
       !
       if (addobs(1, 0, observations(i)%x, observations(i)%y) .lt. 0) then
          gvcart_loaddata_ = -1
          return
       end if
       
    end do
    
    close(FUNIT)

    !
    ! Return 0 for success
    !
    gvcart_loaddata_ = 0

  end function gvcart_loaddata_

  !
  ! This function must compute a predicted value from the model values at the pre-registered
  ! points. That is, for each point, the array values will return the value at that point.
  ! For gradient based sampling (HMC) you need to specify the gradient for the
  ! prediction with respect to each value, ie weights[i] = dprediction/dvalues[i].
  ! For this simple regression problem, prediction = values(1) and hence weight(1) = 1.0.
  !
  ! Special note: the observationidx (not used in the regression case) refers to which observation
  ! this prediction is for in case extra user data is needed for calculation. This index uses
  ! C indexing and hence starts at 0. In fortran a 1 may need to be added to refer to the
  ! correct observation when fortran arrays are used for storage of observations.
  !
  function gvcart_compute_prediction_(nmodels, observationidx, npoints, values, weights, prediction) bind (c)

    integer, intent(in) :: nmodels
    integer, intent(in) :: observationidx
    integer, intent(in) :: npoints
    real(kind = c_double), dimension(npoints), intent(in) :: values
    real(kind = c_double), dimension(npoints), intent(out) :: weights
    real(kind = c_double), intent(out) :: prediction

    integer :: gvcart_compute_prediction_

    weights(1) = 1.0
    prediction = values(1)

    gvcart_compute_prediction_ = 0

  end function gvcart_compute_prediction_

  !
  ! This function computes the likelihood based on the predictions of all observations.
  ! It must compute residuals and weights (d log(like)/dpred[i]). Note that this
  ! computes the negative log likelihood and separates the normalization term
  ! and the exponential term (assuming Gaussian likelihood). The point of this
  ! is that generally the log of the exponential term is more informative without
  ! including the normalization term.
  !
  function gvcart_compute_likelihood_(nmodels, nhierarchical, nobservation, hierarchical, &
       predictions, residuals, weights, like, norm) bind (c)

    integer, intent(in) :: nmodels
    integer, intent(in) :: nhierarchical
    integer, intent(in) :: nobservation

    real(kind = c_double), dimension(nhierarchical), intent(in) :: hierarchical
    real(kind = c_double), dimension(nobservation), intent(in) :: predictions
    real(kind = c_double), dimension(nobservation), intent(out) :: residuals
    real(kind = c_double), dimension(nobservation), intent(out) :: weights
    real(kind = c_double), intent(out) :: like
    real(kind = c_double), intent(out) :: norm

    integer :: gvcart_compute_likelihood_

    !
    ! Local parameters
    !
    real(kind = c_double), parameter :: NORMSCALE = 0.9189385332046727; !0.5*log(2.0*pi)

    real(kind = c_double) :: res
    real(kind = c_double) :: n
    integer :: i
    

    do i = 1, nobservation
       res = predictions(i) - observations(i)%value
       residuals(i) = res
       n = hierarchical(1) * observations(i)%sigma

       weights(i) = res/(n*n);
       
       like = like + (res*weights(i))/(2.0)
       norm = norm + log(n)

    end do

    gvcart_compute_likelihood_ = 0

  end function gvcart_compute_likelihood_

  !
  ! This function is typically used in the generation of synthetic observations
  ! to run synthetic tests. It requires saving of data given by the predictions
  ! array into the format used above in gvcart_loaddata_.
  !
  function gvcart_savedata_(n, filename, noiselevel, nobservation, predictions) bind(c)

    integer, intent(in) :: n
    character, dimension(n), intent(in) :: filename

    real(kind = c_double), intent(in) :: noiselevel
    integer, intent(in) :: nobservation
    real(kind = c_double), dimension(nobservation), intent(in) :: predictions

    integer :: gvcart_savedata_

    !
    ! Local parameters
    !

    integer, parameter :: FUNIT = 10
    
    integer :: status
    integer :: i
    
    open(FUNIT, file = copystring(filename), status = 'replace', action = 'write', iostat = status)
    if (status .ne. 0) then
       gvcart_savedata_ = -1
       return
    end if
    
    write(FUNIT, *) nobservation
    
    do i = 1, nobservation
       
       write(FUNIT, *) observations(i)%x(1), observations(i)%y(1), predictions(i), noiselevel
       
    end do
    
    close(FUNIT)
    
    gvcart_savedata_ = 0
    
  end function gvcart_savedata_

end module genericregression
