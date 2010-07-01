# Microsoft Developer Studio Project File - Name="IFEMfluiddeformable" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=IFEMfluiddeformable - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "IFEMfluiddeformable.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "IFEMfluiddeformable.mak" CFG="IFEMfluiddeformable - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "IFEMfluiddeformable - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "IFEMfluiddeformable - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "IFEMfluiddeformable - Win32 Release"

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

!ELSEIF  "$(CFG)" == "IFEMfluiddeformable - Win32 Debug"

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
# ADD F90 /check:bounds /compile_only /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib /nologo /stack:0x20001ff4,0x10001005 /subsystem:console /incremental:no /debug /machine:I386 /pdbtype:sept

!ENDIF 

# Begin Target

# Name "IFEMfluiddeformable - Win32 Release"
# Name "IFEMfluiddeformable - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Source File

SOURCE=.\block.f90
# End Source File
# Begin Source File

SOURCE=.\blockgmres.f90
# End Source File
# Begin Source File

SOURCE=.\correct.f90
# End Source File
# Begin Source File

SOURCE=.\delta_nonuniform.f90
# End Source File
# Begin Source File

SOURCE=.\echoinput.f90
# End Source File
# Begin Source File

SOURCE=.\ensight_output.f90
# End Source File
# Begin Source File

SOURCE=.\equal.f90
# End Source File
# Begin Source File

SOURCE=.\error.f90
# End Source File
# Begin Source File

SOURCE=.\facemap.f90
# End Source File
# Begin Source File

SOURCE=.\fluid_variables.f90
# End Source File
# Begin Source File

SOURCE=.\form.f90
# End Source File
# Begin Source File

SOURCE=.\gaussj.f90
# End Source File
# Begin Source File

SOURCE=.\gjinv.f90
# End Source File
# Begin Source File

SOURCE=.\global_constants.f90
# End Source File
# Begin Source File

SOURCE=.\global_simulation_parameter.f90
# End Source File
# Begin Source File

SOURCE=.\gmres.f90
# End Source File
# Begin Source File

SOURCE=.\hg.f90
# End Source File
# Begin Source File

SOURCE=.\hypo.f90
# End Source File
# Begin Source File

SOURCE=.\initialize.f90
# End Source File
# Begin Source File

SOURCE=.\length.f90
# End Source File
# Begin Source File

SOURCE=.\main.f90
# End Source File
# Begin Source File

SOURCE=.\meshgen_fluid.f90
# End Source File
# Begin Source File

SOURCE=.\meshgen_solid.f90
# End Source File
# Begin Source File

SOURCE=.\nondimension.f90
# End Source File
# Begin Source File

SOURCE=.\norm.f90
# End Source File
# Begin Source File

SOURCE=.\parseinput.f90
# End Source File
# Begin Source File

SOURCE=.\quad2d3n.f90
# End Source File
# Begin Source File

SOURCE=.\quad2d4n.f90
# End Source File
# Begin Source File

SOURCE=.\quad3d4n.f90
# End Source File
# Begin Source File

SOURCE=.\quad3d8n.f90
# End Source File
# Begin Source File

SOURCE=.\r_bdpd_curr.f90
# End Source File
# Begin Source File

SOURCE=.\r_bdpd_init.f90
# End Source File
# Begin Source File

SOURCE=.\r_common.f90
# End Source File
# Begin Source File

SOURCE=.\r_element.f90
# End Source File
# Begin Source File

SOURCE=.\r_jacob.f90
# End Source File
# Begin Source File

SOURCE=.\r_load.f90
# End Source File
# Begin Source File

SOURCE=.\r_nodalf.f90
# End Source File
# Begin Source File

SOURCE=.\r_sboc.f90
# End Source File
# Begin Source File

SOURCE=.\r_sbpress.f90
# End Source File
# Begin Source File

SOURCE=.\r_scal.f90
# End Source File
# Begin Source File

SOURCE=.\r_scauchy.f90
# End Source File
# Begin Source File

SOURCE=.\r_smaterj.f90
# End Source File
# Begin Source File

SOURCE=.\r_spiola.f90
# End Source File
# Begin Source File

SOURCE=.\r_spiola_elastic.f90
# End Source File
# Begin Source File

SOURCE=.\r_spiola_viscous.f90
# End Source File
# Begin Source File

SOURCE=.\r_spress.f90
# End Source File
# Begin Source File

SOURCE=.\r_sreadinit.f90
# End Source File
# Begin Source File

SOURCE=.\r_sstif.f90
# End Source File
# Begin Source File

SOURCE=.\r_sstrain.f90
# End Source File
# Begin Source File

SOURCE=.\r_stang.f90
# End Source File
# Begin Source File

SOURCE=.\r_stoxc.f90
# End Source File
# Begin Source File

SOURCE=.\r_timefun.f90
# End Source File
# Begin Source File

SOURCE=.\read.f90
# End Source File
# Begin Source File

SOURCE=.\rkpmshape2d.f90
# End Source File
# Begin Source File

SOURCE=.\rkpmshape3d.f90
# End Source File
# Begin Source File

SOURCE=.\run_variables.f90
# End Source File
# Begin Source File

SOURCE=.\set.f90
# End Source File
# Begin Source File

SOURCE=.\shape.f90
# End Source File
# Begin Source File

SOURCE=.\solid_solver.f90
# End Source File
# Begin Source File

SOURCE=.\solid_update.f90
# End Source File
# Begin Source File

SOURCE=.\solid_variables.f90
# End Source File
# Begin Source File

SOURCE=.\update.f90
# End Source File
# Begin Source File

SOURCE=.\vol.f90
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl;fi;fd"
# Begin Source File

SOURCE=.\hypo_declaration_fluid.fi
# End Source File
# Begin Source File

SOURCE=.\hypo_declaration_solid.fi
# End Source File
# Begin Source File

SOURCE=.\hypo_fluid_solver.fi
# End Source File
# Begin Source File

SOURCE=.\hypo_prepare_fluid.fi
# End Source File
# Begin Source File

SOURCE=.\hypo_prepare_solid.fi
# End Source File
# Begin Source File

SOURCE=.\hypo_restart_file_check.fi
# End Source File
# Begin Source File

SOURCE=.\hypo_restart_read.fi
# End Source File
# Begin Source File

SOURCE=.\hypo_restart_write.fi
# End Source File
# Begin Source File

SOURCE=.\hypo_write_output.fi
# End Source File
# Begin Source File

SOURCE=.\vol2d3n.fi
# End Source File
# Begin Source File

SOURCE=.\vol2d4n.fi
# End Source File
# Begin Source File

SOURCE=.\vol3d4n.fi
# End Source File
# Begin Source File

SOURCE=.\vol3d8n.fi
# End Source File
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"
# End Group
# End Target
# End Project
