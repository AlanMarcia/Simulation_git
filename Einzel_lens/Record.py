# ----------------------------------------------
# Script Recorded by Ansys Electronics Desktop Version 2024.2.0
# 15:05:24  Feb 10, 2026
# ----------------------------------------------
import ScriptEnv
ScriptEnv.Initialize("Ansoft.ElectronicsDesktop")
oDesktop.RestoreWindow()
oProject = oDesktop.SetActiveProject("Einzel_lens")
oDesign = oProject.SetActiveDesign("Maxwell2DDesign2")
oModule = oDesign.GetModule("FieldsReporter")
oModule.EnterQty("E")
oModule.ExportOnGrid("D:\\Users\\Public\\Simulation_git\\Einzel_lens\\export_field\\field_neg100.fld", ["-0.38mm", "-0.72mm", "0mm"], ["5.7mm", "0.72mm", "0mm"], ["0.01mm", "0.01mm", "0mm"], "Setup1 : LastAdaptive", 
	[
		"Vbias:="		, "-100V"
	], 
	[
		"NAME:ExportOption",
		"IncludePtInOutput:="	, True,
		"RefCSName:="		, "Global",
		"PtInSI:="		, True,
		"FieldInRefCS:="	, False
	], "Cartesian", ["0mm", "0mm", "0mm"], False)

oModule.ExportOnGrid("D:\\Users\\Public\\Simulation_git\\Einzel_lens\\export_field\\field_neg200.fld", ["-0.38mm", "-0.72mm", "0mm"], ["5.7mm", "0.72mm", "0mm"], ["0.01mm", "0.01mm", "0mm"], "Setup1 : LastAdaptive", 
	[
		"Vbias:="		, "-200V"
	], 
	[
		"NAME:ExportOption",
		"IncludePtInOutput:="	, True,
		"RefCSName:="		, "Global",
		"PtInSI:="		, True,
		"FieldInRefCS:="	, False
	], "Cartesian", ["0mm", "0mm", "0mm"], False)

oModule.ExportOnGrid("D:\\Users\\Public\\Simulation_git\\Einzel_lens\\export_field\\field_neg300.fld", ["-0.38mm", "-0.72mm", "0mm"], ["5.7mm", "0.72mm", "0mm"], ["0.01mm", "0.01mm", "0mm"], "Setup1 : LastAdaptive", 
	[
		"Vbias:="		, "-300V"
	], 
	[
		"NAME:ExportOption",
		"IncludePtInOutput:="	, True,
		"RefCSName:="		, "Global",
		"PtInSI:="		, True,
		"FieldInRefCS:="	, False
	], "Cartesian", ["0mm", "0mm", "0mm"], False)

oModule.ExportOnGrid("D:\\Users\\Public\\Simulation_git\\Einzel_lens\\export_field\\field_neg400.fld", ["-0.38mm", "-0.72mm", "0mm"], ["5.7mm", "0.72mm", "0mm"], ["0.01mm", "0.01mm", "0mm"], "Setup1 : LastAdaptive", 
	[
		"Vbias:="		, "-400V"
	], 
	[
		"NAME:ExportOption",
		"IncludePtInOutput:="	, True,
		"RefCSName:="		, "Global",
		"PtInSI:="		, True,
		"FieldInRefCS:="	, False
	], "Cartesian", ["0mm", "0mm", "0mm"], False)

oModule.ExportOnGrid("D:\\Users\\Public\\Simulation_git\\Einzel_lens\\export_field\\field_neg500.fld", ["-0.38mm", "-0.72mm", "0mm"], ["5.7mm", "0.72mm", "0mm"], ["0.01mm", "0.01mm", "0mm"], "Setup1 : LastAdaptive", 
	[
		"Vbias:="		, "-500V"
	], 
	[
		"NAME:ExportOption",
		"IncludePtInOutput:="	, True,
		"RefCSName:="		, "Global",
		"PtInSI:="		, True,
		"FieldInRefCS:="	, False
	], "Cartesian", ["0mm", "0mm", "0mm"], False)

oModule.ExportOnGrid("D:\\Users\\Public\\Simulation_git\\Einzel_lens\\export_field\\field_neg600.fld", ["-0.38mm", "-0.72mm", "0mm"], ["5.7mm", "0.72mm", "0mm"], ["0.01mm", "0.01mm", "0mm"], "Setup1 : LastAdaptive", 
	[
		"Vbias:="		, "-600V"
	], 
	[
		"NAME:ExportOption",
		"IncludePtInOutput:="	, True,
		"RefCSName:="		, "Global",
		"PtInSI:="		, True,
		"FieldInRefCS:="	, False
	], "Cartesian", ["0mm", "0mm", "0mm"], False)

oModule.ExportOnGrid("D:\\Users\\Public\\Simulation_git\\Einzel_lens\\export_field\\field_neg700.fld", ["-0.38mm", "-0.72mm", "0mm"], ["5.7mm", "0.72mm", "0mm"], ["0.01mm", "0.01mm", "0mm"], "Setup1 : LastAdaptive", 
	[
		"Vbias:="		, "-700V"
	], 
	[
		"NAME:ExportOption",
		"IncludePtInOutput:="	, True,
		"RefCSName:="		, "Global",
		"PtInSI:="		, True,
		"FieldInRefCS:="	, False
	], "Cartesian", ["0mm", "0mm", "0mm"], False)

oModule.ExportOnGrid("D:\\Users\\Public\\Simulation_git\\Einzel_lens\\export_field\\field_neg800.fld", ["-0.38mm", "-0.72mm", "0mm"], ["5.7mm", "0.72mm", "0mm"], ["0.01mm", "0.01mm", "0mm"], "Setup1 : LastAdaptive", 
	[
		"Vbias:="		, "-800V"
	], 
	[
		"NAME:ExportOption",
		"IncludePtInOutput:="	, True,
		"RefCSName:="		, "Global",
		"PtInSI:="		, True,
		"FieldInRefCS:="	, False
	], "Cartesian", ["0mm", "0mm", "0mm"], False)

oModule.ExportOnGrid("D:\\Users\\Public\\Simulation_git\\Einzel_lens\\export_field\\field_neg900.fld", ["-0.38mm", "-0.72mm", "0mm"], ["5.7mm", "0.72mm", "0mm"], ["0.01mm", "0.01mm", "0mm"], "Setup1 : LastAdaptive", 
	[
		"Vbias:="		, "-900V"
	], 
	[
		"NAME:ExportOption",
		"IncludePtInOutput:="	, True,
		"RefCSName:="		, "Global",
		"PtInSI:="		, True,
		"FieldInRefCS:="	, False
	], "Cartesian", ["0mm", "0mm", "0mm"], False)

oModule.ExportOnGrid("D:\\Users\\Public\\Simulation_git\\Einzel_lens\\export_field\\field_neg1000.fld", ["-0.38mm", "-0.72mm", "0mm"], ["5.7mm", "0.72mm", "0mm"], ["0.01mm", "0.01mm", "0mm"], "Setup1 : LastAdaptive", 
	[
		"Vbias:="		, "-1000V"
	], 
	[
		"NAME:ExportOption",
		"IncludePtInOutput:="	, True,
		"RefCSName:="		, "Global",
		"PtInSI:="		, True,
		"FieldInRefCS:="	, False
	], "Cartesian", ["0mm", "0mm", "0mm"], False)

oModule.ExportOnGrid("D:\\Users\\Public\\Simulation_git\\Einzel_lens\\export_field\\field_neg1500.fld", ["-0.38mm", "-0.72mm", "0mm"], ["5.7mm", "0.72mm", "0mm"], ["0.01mm", "0.01mm", "0mm"], "Setup1 : LastAdaptive", 
	[
		"Vbias:="		, "-1500V"
	], 
	[
		"NAME:ExportOption",
		"IncludePtInOutput:="	, True,
		"RefCSName:="		, "Global",
		"PtInSI:="		, True,
		"FieldInRefCS:="	, False
	], "Cartesian", ["0mm", "0mm", "0mm"], False)

oModule.ExportOnGrid("D:\\Users\\Public\\Simulation_git\\Einzel_lens\\export_field\\field_neg2000.fld", ["-0.38mm", "-0.72mm", "0mm"], ["5.7mm", "0.72mm", "0mm"], ["0.01mm", "0.01mm", "0mm"], "Setup1 : LastAdaptive", 
	[
		"Vbias:="		, "-2000V"
	], 
	[
		"NAME:ExportOption",
		"IncludePtInOutput:="	, True,
		"RefCSName:="		, "Global",
		"PtInSI:="		, True,
		"FieldInRefCS:="	, False
	], "Cartesian", ["0mm", "0mm", "0mm"], False)

oModule.ExportOnGrid("D:\\Users\\Public\\Simulation_git\\Einzel_lens\\export_field\\field_neg2500.fld", ["-0.38mm", "-0.72mm", "0mm"], ["5.7mm", "0.72mm", "0mm"], ["0.01mm", "0.01mm", "0mm"], "Setup1 : LastAdaptive", 
	[
		"Vbias:="		, "-2500V"
	], 
	[
		"NAME:ExportOption",
		"IncludePtInOutput:="	, True,
		"RefCSName:="		, "Global",
		"PtInSI:="		, True,
		"FieldInRefCS:="	, False
	], "Cartesian", ["0mm", "0mm", "0mm"], False)

oModule.ExportOnGrid("D:\\Users\\Public\\Simulation_git\\Einzel_lens\\export_field\\field_neg3000.fld", ["-0.38mm", "-0.72mm", "0mm"], ["5.7mm", "0.72mm", "0mm"], ["0.01mm", "0.01mm", "0mm"], "Setup1 : LastAdaptive", 
	[
		"Vbias:="		, "-3000V"
	], 
	[
		"NAME:ExportOption",
		"IncludePtInOutput:="	, True,
		"RefCSName:="		, "Global",
		"PtInSI:="		, True,
		"FieldInRefCS:="	, False
	], "Cartesian", ["0mm", "0mm", "0mm"], False)

oModule.ExportOnGrid("D:\\Users\\Public\\Simulation_git\\Einzel_lens\\export_field\\field_neg3500.fld", ["-0.38mm", "-0.72mm", "0mm"], ["5.7mm", "0.72mm", "0mm"], ["0.01mm", "0.01mm", "0mm"], "Setup1 : LastAdaptive", 
	[
		"Vbias:="		, "-3500V"
	], 
	[
		"NAME:ExportOption",
		"IncludePtInOutput:="	, True,
		"RefCSName:="		, "Global",
		"PtInSI:="		, True,
		"FieldInRefCS:="	, False
	], "Cartesian", ["0mm", "0mm", "0mm"], False)

oModule.ExportOnGrid("D:\\Users\\Public\\Simulation_git\\Einzel_lens\\export_field\\field_neg4000.fld", ["-0.38mm", "-0.72mm", "0mm"], ["5.7mm", "0.72mm", "0mm"], ["0.01mm", "0.01mm", "0mm"], "Setup1 : LastAdaptive", 
	[
		"Vbias:="		, "-4000V"
	], 
	[
		"NAME:ExportOption",
		"IncludePtInOutput:="	, True,
		"RefCSName:="		, "Global",
		"PtInSI:="		, True,
		"FieldInRefCS:="	, False
	], "Cartesian", ["0mm", "0mm", "0mm"], False)

oModule.ExportOnGrid("D:\\Users\\Public\\Simulation_git\\Einzel_lens\\export_field\\field_neg4500.fld", ["-0.38mm", "-0.72mm", "0mm"], ["5.7mm", "0.72mm", "0mm"], ["0.01mm", "0.01mm", "0mm"], "Setup1 : LastAdaptive", 
	[
		"Vbias:="		, "-4500V"
	], 
	[
		"NAME:ExportOption",
		"IncludePtInOutput:="	, True,
		"RefCSName:="		, "Global",
		"PtInSI:="		, True,
		"FieldInRefCS:="	, False
	], "Cartesian", ["0mm", "0mm", "0mm"], False)

oModule.ExportOnGrid("D:\\Users\\Public\\Simulation_git\\Einzel_lens\\export_field\\field_neg5000.fld", ["-0.38mm", "-0.72mm", "0mm"], ["5.7mm", "0.72mm", "0mm"], ["0.01mm", "0.01mm", "0mm"], "Setup1 : LastAdaptive", 
	[
		"Vbias:="		, "-5000V"
	], 
	[
		"NAME:ExportOption",
		"IncludePtInOutput:="	, True,
		"RefCSName:="		, "Global",
		"PtInSI:="		, True,
		"FieldInRefCS:="	, False
	], "Cartesian", ["0mm", "0mm", "0mm"], False)

oModule.ExportOnGrid("D:\\Users\\Public\\Simulation_git\\Einzel_lens\\export_field\\field_neg5500.fld", ["-0.38mm", "-0.72mm", "0mm"], ["5.7mm", "0.72mm", "0mm"], ["0.01mm", "0.01mm", "0mm"], "Setup1 : LastAdaptive", 
	[
		"Vbias:="		, "-5500V"
	], 
	[
		"NAME:ExportOption",
		"IncludePtInOutput:="	, True,
		"RefCSName:="		, "Global",
		"PtInSI:="		, True,
		"FieldInRefCS:="	, False
	], "Cartesian", ["0mm", "0mm", "0mm"], False)

oModule.ExportOnGrid("D:\\Users\\Public\\Simulation_git\\Einzel_lens\\export_field\\field_neg6000.fld", ["-0.38mm", "-0.72mm", "0mm"], ["5.7mm", "0.72mm", "0mm"], ["0.01mm", "0.01mm", "0mm"], "Setup1 : LastAdaptive", 
	[
		"Vbias:="		, "-6000V"
	], 
	[
		"NAME:ExportOption",
		"IncludePtInOutput:="	, True,
		"RefCSName:="		, "Global",
		"PtInSI:="		, True,
		"FieldInRefCS:="	, False
	], "Cartesian", ["0mm", "0mm", "0mm"], False)

oModule.ExportOnGrid("D:\\Users\\Public\\Simulation_git\\Einzel_lens\\export_field\\field_neg6500.fld", ["-0.38mm", "-0.72mm", "0mm"], ["5.7mm", "0.72mm", "0mm"], ["0.01mm", "0.01mm", "0mm"], "Setup1 : LastAdaptive", 
	[
		"Vbias:="		, "-6500V"
	], 
	[
		"NAME:ExportOption",
		"IncludePtInOutput:="	, True,
		"RefCSName:="		, "Global",
		"PtInSI:="		, True,
		"FieldInRefCS:="	, False
	], "Cartesian", ["0mm", "0mm", "0mm"], False)

oModule.ExportOnGrid("D:\\Users\\Public\\Simulation_git\\Einzel_lens\\export_field\\field_neg7000.fld", ["-0.38mm", "-0.72mm", "0mm"], ["5.7mm", "0.72mm", "0mm"], ["0.01mm", "0.01mm", "0mm"], "Setup1 : LastAdaptive", 
	[
		"Vbias:="		, "-7000V"
	], 
	[
		"NAME:ExportOption",
		"IncludePtInOutput:="	, True,
		"RefCSName:="		, "Global",
		"PtInSI:="		, True,
		"FieldInRefCS:="	, False
	], "Cartesian", ["0mm", "0mm", "0mm"], False)

oModule.ExportOnGrid("D:\\Users\\Public\\Simulation_git\\Einzel_lens\\export_field\\field_neg7500.fld", ["-0.38mm", "-0.72mm", "0mm"], ["5.7mm", "0.72mm", "0mm"], ["0.01mm", "0.01mm", "0mm"], "Setup1 : LastAdaptive", 
	[
		"Vbias:="		, "-7500V"
	], 
	[
		"NAME:ExportOption",
		"IncludePtInOutput:="	, True,
		"RefCSName:="		, "Global",
		"PtInSI:="		, True,
		"FieldInRefCS:="	, False
	], "Cartesian", ["0mm", "0mm", "0mm"], False)

oModule.ExportOnGrid("D:\\Users\\Public\\Simulation_git\\Einzel_lens\\export_field\\field_neg8000.fld", ["-0.38mm", "-0.72mm", "0mm"], ["5.7mm", "0.72mm", "0mm"], ["0.01mm", "0.01mm", "0mm"], "Setup1 : LastAdaptive", 
	[
		"Vbias:="		, "-8000V"
	], 
	[
		"NAME:ExportOption",
		"IncludePtInOutput:="	, True,
		"RefCSName:="		, "Global",
		"PtInSI:="		, True,
		"FieldInRefCS:="	, False
	], "Cartesian", ["0mm", "0mm", "0mm"], False)

oModule.ExportOnGrid("D:\\Users\\Public\\Simulation_git\\Einzel_lens\\export_field\\field_neg8500.fld", ["-0.38mm", "-0.72mm", "0mm"], ["5.7mm", "0.72mm", "0mm"], ["0.01mm", "0.01mm", "0mm"], "Setup1 : LastAdaptive", 
	[
		"Vbias:="		, "-8500V"
	], 
	[
		"NAME:ExportOption",
		"IncludePtInOutput:="	, True,
		"RefCSName:="		, "Global",
		"PtInSI:="		, True,
		"FieldInRefCS:="	, False
	], "Cartesian", ["0mm", "0mm", "0mm"], False)

oModule.ExportOnGrid("D:\\Users\\Public\\Simulation_git\\Einzel_lens\\export_field\\field_neg9000.fld", ["-0.38mm", "-0.72mm", "0mm"], ["5.7mm", "0.72mm", "0mm"], ["0.01mm", "0.01mm", "0mm"], "Setup1 : LastAdaptive", 
	[
		"Vbias:="		, "-9000V"
	], 
	[
		"NAME:ExportOption",
		"IncludePtInOutput:="	, True,
		"RefCSName:="		, "Global",
		"PtInSI:="		, True,
		"FieldInRefCS:="	, False
	], "Cartesian", ["0mm", "0mm", "0mm"], False)

oModule.ExportOnGrid("D:\\Users\\Public\\Simulation_git\\Einzel_lens\\export_field\\field_neg9500.fld", ["-0.38mm", "-0.72mm", "0mm"], ["5.7mm", "0.72mm", "0mm"], ["0.01mm", "0.01mm", "0mm"], "Setup1 : LastAdaptive", 
	[
		"Vbias:="		, "-9500V"
	], 
	[
		"NAME:ExportOption",
		"IncludePtInOutput:="	, True,
		"RefCSName:="		, "Global",
		"PtInSI:="		, True,
		"FieldInRefCS:="	, False
	], "Cartesian", ["0mm", "0mm", "0mm"], False)

oModule.ExportOnGrid("D:\\Users\\Public\\Simulation_git\\Einzel_lens\\export_field\\field_neg10000.fld", ["-0.38mm", "-0.72mm", "0mm"], ["5.7mm", "0.72mm", "0mm"], ["0.01mm", "0.01mm", "0mm"], "Setup1 : LastAdaptive", 
	[
		"Vbias:="		, "-10000V"
	], 
	[
		"NAME:ExportOption",
		"IncludePtInOutput:="	, True,
		"RefCSName:="		, "Global",
		"PtInSI:="		, True,
		"FieldInRefCS:="	, False
	], "Cartesian", ["0mm", "0mm", "0mm"], False)