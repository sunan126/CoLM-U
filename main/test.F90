#include <define.h>

  SUBROUTINE LuLccEnergyConserve ()

    ! 需要什么变量的话就要USE那个定义了该变量的F90文件
    USE MOD_LuLccTimeInvars
    USE MOD_LuLccTimeVars
    USE MOD_TimeVariables

    IMPLICIT NONE
    
    INTEGER, allocatable :: tmp(:) !某个patch的lc类型在to序列中的位置
    INTEGER, allocatable :: tmpfr(:) !某个patch的lc来源
    INTEGER, allocatable :: lccpct(:) !某个patch的lc来源百分比，包括自己
    INTEGER, allocatable :: tmpnp_(:) !某个patch的来源patch的索引
    INTEGER :: fr(nnum) !lc旧序列
    INTEGER :: to(nnum) !lc新序列
    INTEGER :: K

     	DO j = 1, lat_points
        DO i = 1, lon_points
          np = grid_patch_s (i,j) !初始化后的编号
          np_= grid_patch_s_(i,j) !保存了上一年的编号


#ifdef IGBP_CLASSIFICATION

          fr = lcfr(:,np)
          to = lcto(:,np)

          DO WHILE (np.le.grid_patch_e(i,j))
          ! 编号不能超过该经纬度网格对应的一串编号
            num = count(to. eq. patchclass(np))
            ! 统计: to序列中lc==np的个数,并获取索引号 ->继而获取from序列的lc、LCC百分比
            IF (num .ge. 1) THEN
              allocate (tmp(num))
              allocate (tmpfr(num))
              allocate (lccpct(num))

              tmp(1:1) = findloc(to, patchclass(np))
              tmpfr(1) = fr(tmp(1))
              lccpct(1) = lccpct(tmp(1))
                              
              DO k = 1, num
                
                IF (k .ge. 2) THEN
                  tmp(k:k) = findloc(to(tmp(k-1)+1:50),patchclass(np)) + tmp(k-1)
                  tmpfr(k) = fr(tmp(k))
                  lccpct(k) = lccpct(tmp(k))
                ENDIF
                
                DO WHILE (np_ .le. grid_patch_e_(i,j))
                  ! 根据 tmpfr获取该类型的索引，然后存到一个动态数组
                  IF (patchclass(np_) .eq. tmpfr(k)) THEN
                    tmpnp_(k) = np_
                  ENDIF
                ENDDO

              ENDDO

              ! 计算要加权平均的变量
              wliq_soisno (:,np) = 0
              DO k = 1, num
                wliq_soisno (:,np) = wliq_soisno (:,np) + wliq_soisno_(:,tmpnp_(k))*lccpct(k))/sum(lccpct(1:num))
              ENDDO

            ELSE
            ! 说明新一年的这种类型并不包含在50对lcc的汇中，即这这种地表占比非常小。
              print*, 'patch with little fraction'

              DO WHILE (np_ .le. grid_patch_e_(i,j))
                ! 如果去年存在同样的地表，直接赋值；如果去年不存在同样的地表，因为占比很小，就直接用初始化的值
                IF (patchclass(np_) .eq. patchclass(np)) THEN
                  wliq_soisno (:,np) = wliq_soisno_(:,np_)
                  print*, 'old patch', 这种地表的占比.
                  EXIT
                ENDIF
              ENDDO

            ENDIF
            
            np = np + 1
          ENDDO
            

#endif

END SUBROUTINE LuLccEnergyConserve