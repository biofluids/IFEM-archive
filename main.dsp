# Microsoft Developer Studio Project File - Name="main" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=main - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "main.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "main.mak" CFG="main - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "main - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "main - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "main - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE F90 /compile_only /nologo /warn:nofileopt
# ADD F90 /compile_only /nologo /warn:nofileopt
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib /nologo /subsystem:console /machine:I386

!ELSEIF  "$(CFG)" == "main - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE F90 /check:bounds /compile_only /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD F90 /browser /check:bounds /compile_only /debug:full /nologo /optimize:5 /threads /traceback /warn:argument_checking /warn:nofileopt
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /ML /W4 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /FR /YX /FD /GZ /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept

!ENDIF 

# Begin Target

# Name "main - Win32 Release"
# Name "main - Win32 Debug"
# Begin Source File

SOURCE=.\block.f
DEP_F90_BLOCK=\
	".\global.h"\
	".\sh3d4n.h"\
	".\sh3d8n.h"\
	

!IF  "$(CFG)" == "main - Win32 Release"

!ELSEIF  "$(CFG)" == "main - Win32 Debug"

# ADD F90 /threads

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\blockgmres.f
DEP_F90_BLOCKG=\
	".\global.h"\
	".\sh3d4n.h"\
	".\sh3d8n.h"\
	
# End Source File
# Begin Source File

SOURCE=.\calcaccel.f
DEP_F90_CALCA=\
	".\main_common"\
	
# End Source File
# Begin Source File

SOURCE=.\correct3dl.f
# End Source File
# Begin Source File

SOURCE=.\declaration_fluid.fi
# End Source File
# Begin Source File

SOURCE=.\declaration_solid.fi
# End Source File
# Begin Source File

SOURCE=.\delta.f
# End Source File
# Begin Source File

SOURCE=.\delta_nonuniform.f
DEP_F90_DELTA=\
	".\global.h"\
	".\vol3d4n.h"\
	".\vol3d8n.h"\
	
# End Source File
# Begin Source File

SOURCE=.\disturbance.f
DEP_F90_DISTU=\
	".\main_common"\
	".\r_common"\
	
# End Source File
# Begin Source File

SOURCE=.\echoinput.f
DEP_F90_ECHOI=\
	".\global.h"\
	
# End Source File
# Begin Source File

SOURCE=.\equal.f
# End Source File
# Begin Source File

SOURCE=.\error.f
# End Source File
# Begin Source File

SOURCE=.\f_fiber1.f
# End Source File
# Begin Source File

SOURCE=.\f_fiber2.f
# End Source File
# Begin Source File

SOURCE=.\facemap.f
DEP_F90_FACEM=\
	".\global.h"\
	
# End Source File
# Begin Source File

SOURCE=.\fclear.f
# End Source File
# Begin Source File

SOURCE=.\fem_fluid_solver.fi
# End Source File
# Begin Source File

SOURCE=.\fiber1.f
DEP_F90_FIBER=\
	".\main_common"\
	
# End Source File
# Begin Source File

SOURCE=.\fiber10.f
DEP_F90_FIBER1=\
	".\main_common"\
	
# End Source File
# Begin Source File

SOURCE=.\fiber11.f
DEP_F90_FIBER11=\
	".\main_common"\
	".\r_common"\
	
# End Source File
# Begin Source File

SOURCE=.\fiber13.f
DEP_F90_FIBER13=\
	".\iba_application_parameters.fh"\
	".\iba_application_variables.fh"\
	".\main_common"\
	
# End Source File
# Begin Source File

SOURCE=.\fiber2.f
DEP_F90_FIBER2=\
	".\iba_application_parameters.fh"\
	".\iba_application_variables.fh"\
	
# End Source File
# Begin Source File

SOURCE=.\fiber3.f
DEP_F90_FIBER3=\
	".\main_common"\
	".\r_common"\
	
# End Source File
# Begin Source File

SOURCE=.\fiber5.f
DEP_F90_FIBER5=\
	".\main_common"\
	
# End Source File
# Begin Source File

SOURCE=.\fiber8.f
DEP_F90_FIBER8=\
	".\main_common"\
	
# End Source File
# Begin Source File

SOURCE=.\fiber9.f
DEP_F90_FIBER9=\
	".\main_common"\
	
# End Source File
# Begin Source File

SOURCE=.\form.f
DEP_F90_FORM_=\
	".\global.h"\
	".\malloc.h"\
	
# End Source File
# Begin Source File

SOURCE=.\gaussj.f
# End Source File
# Begin Source File

SOURCE=.\gjinv.f
# End Source File
# Begin Source File

SOURCE=.\gmres.f
DEP_F90_GMRES=\
	".\global.h"\
	
# End Source File
# Begin Source File

SOURCE=.\hg.f
DEP_F90_HG_F24=\
	".\global.h"\
	
# End Source File
# Begin Source File

SOURCE=.\hydro.f
DEP_F90_HYDRO=\
	".\global.h"\
	
# End Source File
# Begin Source File

SOURCE=.\hypo.f
DEP_F90_HYPO_=\
	".\declaration_fluid.fi"\
	".\declaration_solid.fi"\
	".\fem_fluid_solver.fi"\
	".\global.h"\
	".\iba_application_parameters.fh"\
	".\ibd0_app_vals_nonzerocfddiv.fh"\
	".\ibd0_exchange_pars.fh"\
	".\ibd0_nonzerocfddiv_pars.fh"\
	".\ibd0_nonzerocfddiv_vars.fh"\
	".\ibg_change_me_ptcon_var_common_equiv.fh"\
	".\ibg_change_me_ptcon_var_decl.fh"\
	".\ibg_parameters_run.fh"\
	".\ibg_variables_cfd.fh"\
	".\ibg_variables_cloud.fh"\
	".\ibg_variables_domain.fh"\
	".\ibg_variables_point.fh"\
	".\ibg_variables_run.fh"\
	".\main_common"\
	".\malloc.h"\
	".\pointer.fi"\
	".\prepare_fluid.fi"\
	".\prepare_solid.fi"\
	".\r_common"\
	".\solids_solver.fi"\
	".\solids_update.fi"\
	".\write_output.fi"\
	
# End Source File
# Begin Source File

SOURCE=.\initialize.f
DEP_F90_INITI=\
	".\global.h"\
	
# End Source File
# Begin Source File

SOURCE=.\io8.f
DEP_F90_IO8_F=\
	".\main_common"\
	".\r_common"\
	
# End Source File
# Begin Source File

SOURCE=.\lenght.f
DEP_F90_LENGH=\
	".\global.h"\
	".\sh3d4n.h"\
	".\sh3d8n.h"\
	
# End Source File
# Begin Source File

SOURCE=.\link2.f
# End Source File
# Begin Source File

SOURCE=.\link3.f
# End Source File
# Begin Source File

SOURCE=.\link4.f
# End Source File
# Begin Source File

SOURCE=.\link6.f
# End Source File
# Begin Source File

SOURCE=.\main.f
DEP_F90_MAIN_=\
	".\global.h"\
	".\main_common"\
	".\r_common"\
	
# End Source File
# Begin Source File

SOURCE=.\meshgen.f
DEP_F90_MESHG=\
	".\global.h"\
	
# End Source File
# Begin Source File

SOURCE=.\movepoints.f
DEP_F90_MOVEP=\
	".\main_common"\
	".\r_common"\
	
# End Source File
# Begin Source File

SOURCE=.\nondimension.f
DEP_F90_NONDI=\
	".\global.h"\
	
# End Source File
# Begin Source File

SOURCE=.\norm.f
DEP_F90_NORM_=\
	".\global.h"\
	
# End Source File
# Begin Source File

SOURCE=.\parseinput_fluid.f
DEP_F90_PARSE=\
	".\global.h"\
	
# End Source File
# Begin Source File

SOURCE=.\pointer.fi
# End Source File
# Begin Source File

SOURCE=.\prepare_fluid.fi
# End Source File
# Begin Source File

SOURCE=.\prepare_solid.fi
# End Source File
# Begin Source File

SOURCE=.\quad2d3n.f
# End Source File
# Begin Source File

SOURCE=.\quad2d4n.f
# End Source File
# Begin Source File

SOURCE=.\quad3d4n.f
# End Source File
# Begin Source File

SOURCE=.\quad3d8n.f
# End Source File
# Begin Source File

SOURCE=.\r_bdpd.f
DEP_F90_R_BDP=\
	".\r_common"\
	
# End Source File
# Begin Source File

SOURCE=.\r_element.f
DEP_F90_R_ELE=\
	".\r_common"\
	
# End Source File
# Begin Source File

SOURCE=.\r_input.f
DEP_F90_R_INP=\
	".\iba_application_parameters.fh"\
	".\ibd0_exchange_pars.fh"\
	".\ibg_change_me_ptcon_var_common_equiv.fh"\
	".\ibg_change_me_ptcon_var_decl.fh"\
	".\ibg_parameters_run.fh"\
	".\ibg_variable_equivalences.fh"\
	".\ibg_variables_cfd.fh"\
	".\ibg_variables_cloud.fh"\
	".\ibg_variables_domain.fh"\
	".\ibg_variables_fluid.fh"\
	".\ibg_variables_marker.fh"\
	".\ibg_variables_point.fh"\
	".\ibg_variables_run.fh"\
	".\main_common"\
	".\r_common"\
	
# End Source File
# Begin Source File

SOURCE=.\r_jacob.f
DEP_F90_R_JAC=\
	".\r_common"\
	
# End Source File
# Begin Source File

SOURCE=.\r_load.f
DEP_F90_R_LOA=\
	".\r_common"\
	
# End Source File
# Begin Source File

SOURCE=.\r_main.f
DEP_F90_R_MAI=\
	".\iba_application_parameters.fh"\
	".\ibd0_exchange_pars.fh"\
	".\ibg_change_me_ptcon_var_common_equiv.fh"\
	".\ibg_change_me_ptcon_var_decl.fh"\
	".\ibg_parameters_run.fh"\
	".\ibg_variable_equivalences.fh"\
	".\ibg_variables_cfd.fh"\
	".\ibg_variables_cloud.fh"\
	".\ibg_variables_domain.fh"\
	".\ibg_variables_fluid.fh"\
	".\ibg_variables_marker.fh"\
	".\ibg_variables_point.fh"\
	".\main_common"\
	".\r_common"\
	
# End Source File
# Begin Source File

SOURCE=.\r_main2.f
DEP_F90_R_MAIN=\
	".\main_common"\
	".\r_common"\
	
# End Source File
# Begin Source File

SOURCE=.\r_nodalf.f
DEP_F90_R_NOD=\
	".\r_common"\
	
# End Source File
# Begin Source File

SOURCE=.\r_print.f
DEP_F90_R_PRI=\
	".\main_common"\
	".\r_common"\
	
# End Source File
# Begin Source File

SOURCE=.\r_sboc.f
DEP_F90_R_SBO=\
	".\r_common"\
	
# End Source File
# Begin Source File

SOURCE=.\r_sbpress.f
DEP_F90_R_SBP=\
	".\r_common"\
	
# End Source File
# Begin Source File

SOURCE=.\r_scalfp.f
DEP_F90_R_SCA=\
	".\r_common"\
	
# End Source File
# Begin Source File

SOURCE=.\r_scalfu.f
DEP_F90_R_SCAL=\
	".\r_common"\
	
# End Source File
# Begin Source File

SOURCE=.\r_scalkpp.f
DEP_F90_R_SCALK=\
	".\r_common"\
	
# End Source File
# Begin Source File

SOURCE=.\r_scalkup.f
DEP_F90_R_SCALKU=\
	".\r_common"\
	
# End Source File
# Begin Source File

SOURCE=.\r_scauchy.f
DEP_F90_R_SCAU=\
	".\r_common"\
	
# End Source File
# Begin Source File

SOURCE=.\r_sinit.f
DEP_F90_R_SIN=\
	".\r_common"\
	
# End Source File
# Begin Source File

SOURCE=.\r_smaterj.f
DEP_F90_R_SMA=\
	".\r_common"\
	
# End Source File
# Begin Source File

SOURCE=.\r_spiola.f
DEP_F90_R_SPI=\
	".\r_common"\
	
# End Source File
# Begin Source File

SOURCE=.\r_spress.f
DEP_F90_R_SPR=\
	".\r_common"\
	
# End Source File
# Begin Source File

SOURCE=.\r_sreadinit.f
DEP_F90_R_SRE=\
	".\r_common"\
	
# End Source File
# Begin Source File

SOURCE=.\r_sstif.f
DEP_F90_R_SST=\
	".\main_common"\
	".\r_common"\
	
# End Source File
# Begin Source File

SOURCE=.\r_sstrain.f
DEP_F90_R_SSTR=\
	".\r_common"\
	
# End Source File
# Begin Source File

SOURCE=.\r_stang.f
DEP_F90_R_STA=\
	".\r_common"\
	
# End Source File
# Begin Source File

SOURCE=.\r_stoxc.f
# End Source File
# Begin Source File

SOURCE=.\r_timefun.f
DEP_F90_R_TIM=\
	".\r_common"\
	
# End Source File
# Begin Source File

SOURCE=.\rkpmshape3d.f
# End Source File
# Begin Source File

SOURCE=.\set.f
DEP_F90_SET_F=\
	".\global.h"\
	
# End Source File
# Begin Source File

SOURCE=.\shape.f
DEP_F90_SHAPE=\
	".\global.h"\
	
# End Source File
# Begin Source File

SOURCE=.\sharp.f
DEP_F90_SHARP=\
	".\global.h"\
	
# End Source File
# Begin Source File

SOURCE=.\solids_solver.fi
# End Source File
# Begin Source File

SOURCE=.\solids_update.fi
# End Source File
# Begin Source File

SOURCE=.\update.f
DEP_F90_UPDAT=\
	".\global.h"\
	
# End Source File
# Begin Source File

SOURCE=.\vol.f
DEP_F90_VOL_F=\
	".\global.h"\
	".\sh3d4n.h"\
	".\sh3d8n.h"\
	

!IF  "$(CFG)" == "main - Win32 Release"

!ELSEIF  "$(CFG)" == "main - Win32 Debug"

# ADD F90 /optimize:0

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\write_output.fi
# End Source File
# Begin Source File

SOURCE=.\zfem_ensCase.f
DEP_F90_ZFEM_=\
	".\main_common"\
	
# End Source File
# Begin Source File

SOURCE=.\zfem_ensFluid.f
DEP_F90_ZFEM_E=\
	".\global.h"\
	".\main_common"\
	".\r_common"\
	
# End Source File
# Begin Source File

SOURCE=.\zfem_ensGeo.f
DEP_F90_ZFEM_EN=\
	".\global.h"\
	".\main_common"\
	".\r_common"\
	
# End Source File
# Begin Source File

SOURCE=.\zfem_tec.f
DEP_F90_ZFEM_T=\
	".\global.h"\
	".\main_common"\
	".\r_common"\
	
# End Source File
# Begin Source File

SOURCE=.\zibm_ensCase.f
DEP_F90_ZIBM_=\
	".\main_common"\
	
# End Source File
# Begin Source File

SOURCE=.\zibm_ensFluid.f
DEP_F90_ZIBM_E=\
	".\global.h"\
	
# End Source File
# Begin Source File

SOURCE=.\zibm_ensGeo.f
DEP_F90_ZIBM_EN=\
	".\global.h"\
	".\main_common"\
	
# End Source File
# Begin Source File

SOURCE=.\zibm_tec.f
DEP_F90_ZIBM_T=\
	".\global.h"\
	".\main_common"\
	
# End Source File
# End Target
# End Project
