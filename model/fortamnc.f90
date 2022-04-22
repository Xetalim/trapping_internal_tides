PROGRAM fortamnc
  USE Fields
  IMPLICIT NONE
  CHARACTER(len=*), PARAMETER :: filename = "/nobackup_1/users/drijfhou/iwa/iwa044/iwa044-12-harmonie"
  INTEGER,          PARAMETER :: id=53,ii=321,jd=322,jj=52,kd=42,kk=44,k1=43

  REAL(kind=4) :: data(id,jd,kd)
  REAL(kind=8) :: datb(ii,jj,kk)

  CALL readInputFile
  CALL rearrangeData
  CALL writeOutputFile
  STOP
CONTAINS

  SUBROUTINE readInputFile
    INTEGER, PARAMETER :: unit = 99
    OPEN(unit=unit,file=filename,form='unformatted',access='sequential')
    READ(unit=unit) data(:,:,:)
    CLOSE(unit=unit)
    RETURN
  END SUBROUTINE readInputFile

  SUBROUTINE rearrangeData
    REAL(kind=4) :: tmp(id,jd,kd)
    INTEGER :: i, j, k
    tmp = data
     DO i=1,id-1
     DO j=1,jd-1
     datb(j,i,1)=0.
     datb(j,i,kk)=0.
     END DO
     END DO
    DO k=1,kd
    DO i=1,id-1
       DO j=1,jd-1
    datb(j,id-i,k+1)=tmp(i,j,k)
       END DO
    END DO
    END DO
    RETURN
  END SUBROUTINE rearrangeData

  SUBROUTINE writeOutputFile
    TYPE(DataFile)   :: file
    TYPE(Error)      :: err
    TYPE(Coordinate) :: xcoord, ycoord, zcoord
    TYPE(Grid3D)     :: xyzgrid
    TYPE(Field3D)    :: data_field
    INTEGER          :: i, j, k

    CALL initialize(file,name="amp_ex44_12.nc",overwrite=.TRUE.)
    CALL initialize(xcoord,name="xc",description="Distance",unit="kilometer",values=(/ (0.0d0+(1200.*(i-1))/(ii-1),i=1,ii) /))
    CALL initialize(ycoord,name="yc",description="Distance",unit="kilometer",values=(/ (0.0d0+(191.25*(j-1))/(jj-1),j=1,jj) /))
    CALL initialize(zcoord,name="zc",description="Depth",unit="meter",values=(/ (0d0- (4300.0*(k-1))/(kk-1),k=1,kk) /))
    CALL initialize(xyzgrid,length=xcoord,width=ycoord,height=zcoord)
    CALL initialize(data_field,name="amplitude",description="Amplitude",unit="",grid=xyzgrid,values=datb(:,:,:))
    CALL write(data_field,file)
    CALL finalize(data_field)
    CALL finalize(xyzgrid)
    CALL finalize(zcoord)
    CALL finalize(ycoord)
    CALL finalize(xcoord)
    CALL finalize(file)
    RETURN
  END SUBROUTINE writeOutputFile

END PROGRAM fortamnc
