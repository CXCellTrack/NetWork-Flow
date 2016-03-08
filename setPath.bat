@echo 添加 ssvm train test 用户环境变量
@echo off
rem 注意并非用户变量！
rem 注意xxx=XXX中间一定不能加空格！否则会报错！

rem set regpath=HKEY_LOCAL_MACHINE\SYSTEM\CurrentControlSet\Control\Session Manager\Environment
set regpath=HKEY_CURRENT_USER\Environment

set evname=trainSegPath

set trainSegPath=E:\datasets\first_edition\training_datasets\N2DL-HeLa\01_2-16_seg

reg add "%regpath%" /v %evname% /d %trainSegPath% /f

set evname=trainTrackPath

set trainTrackPath=E:\datasets\first_edition\training_datasets\N2DL-HeLa\01_2-16_track

reg add "%regpath%" /v %evname% /d %trainTrackPath% /f



set evname=testSegPath

set testSegPath=E:\datasets\first_edition\training_datasets\N2DL-HeLa\02_4-16_seg

reg add "%regpath%" /v %evname% /d %testSegPath% /f

set evname=testTrackPath

set testTrackPath=E:\datasets\first_edition\training_datasets\N2DL-HeLa\02_4-16_track

reg add "%regpath%" /v %evname% /d %testTrackPath% /f






rem 删除系统变量 reg delete "%regpath%" /v "%evname%"  /f

rem 重启任务管理器，用于刷新注册表

rem taskkill /f /im explorer.exe && explorer.exe 
 
rem 暂停但不显示“请按任意键继续”

pause>nul


