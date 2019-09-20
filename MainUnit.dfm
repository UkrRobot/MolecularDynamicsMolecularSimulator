object MainForm: TMainForm
  Left = 238
  Top = 122
  BorderStyle = bsSingle
  Caption = 'Molecular Dynamics Molecular Simulator'
  ClientHeight = 341
  ClientWidth = 690
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'MS Sans Serif'
  Font.Style = []
  KeyPreview = True
  OldCreateOrder = False
  Position = poDesktopCenter
  OnCreate = FormCreate
  OnDestroy = FormDestroy
  PixelsPerInch = 96
  TextHeight = 13
  object ImagePanel: TPanel
    Left = 10
    Top = 10
    Width = 321
    Height = 321
    Hint = 'These are the moleculs'
    BevelOuter = bvLowered
    Color = clWhite
    ParentShowHint = False
    ShowHint = True
    TabOrder = 0
    object AtomView: TImage
      Left = 0
      Top = 0
      Width = 321
      Height = 321
      ParentShowHint = False
      ShowHint = False
    end
  end
  object StartButton: TButton
    Left = 340
    Top = 300
    Width = 71
    Height = 31
    Hint = 'Start/Stop simulation'
    Caption = 'Go'
    ParentShowHint = False
    ShowHint = True
    TabOrder = 1
    OnClick = StartButtonClick
  end
  object ResetButton: TButton
    Left = 500
    Top = 300
    Width = 71
    Height = 31
    Hint = 'Reset the position of atoms'
    Caption = 'Reset'
    ParentShowHint = False
    ShowHint = True
    TabOrder = 2
    OnClick = ResetButtonClick
  end
  object QuitButton: TButton
    Left = 610
    Top = 300
    Width = 71
    Height = 31
    Hint = 'Close this application'
    Cancel = True
    Caption = 'Quit'
    ModalResult = 2
    ParentShowHint = False
    ShowHint = True
    TabOrder = 3
    OnClick = QuitButtonClick
  end
  object OneFrameButton: TButton
    Left = 420
    Top = 300
    Width = 71
    Height = 31
    Hint = 'Simulate one frame'
    Caption = 'One Frame'
    ParentShowHint = False
    ShowHint = True
    TabOrder = 4
    OnClick = OneFrameButtonClick
  end
  object AdditionalPages: TPageControl
    Left = 340
    Top = 10
    Width = 340
    Height = 281
    ActivePage = ParametersTab
    DragKind = dkDock
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -11
    Font.Name = 'MS Sans Serif'
    Font.Style = []
    ParentFont = False
    TabIndex = 0
    TabOrder = 5
    OnChange = AdditionalPagesChange
    object ParametersTab: TTabSheet
      Caption = 'Parameters'
      object PeriodicBorderSwitch: TSpeedButton
        Left = 230
        Top = 220
        Width = 91
        Height = 31
        Hint = 'Set periodic boundary conditions'
        GroupIndex = 1
        Caption = 'Periodic'
        ParentShowHint = False
        ShowHint = True
        OnClick = PeriodicBorderSwitchClick
      end
      object FixedBorderSwitch: TSpeedButton
        Left = 130
        Top = 220
        Width = 91
        Height = 31
        Hint = 'Set fixed boundary conditions'
        GroupIndex = 1
        Down = True
        Caption = 'Fixed'
        ParentShowHint = False
        ShowHint = True
        OnClick = FixedBorderSwitchClick
      end
      object Potential_ELabel: TPanel
        Left = 0
        Top = 40
        Width = 121
        Height = 21
        BevelOuter = bvNone
        Caption = 'Potential depth, K'
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -11
        Font.Name = 'MS Sans Serif'
        Font.Style = [fsBold]
        ParentFont = False
        TabOrder = 0
      end
      object Potential_DLabel: TPanel
        Left = 0
        Top = 70
        Width = 121
        Height = 21
        BevelOuter = bvNone
        Caption = 'Potential minimum'
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -11
        Font.Name = 'MS Sans Serif'
        Font.Style = [fsBold]
        ParentFont = False
        TabOrder = 1
      end
      object Potential_EEdit: TPanel
        Left = 130
        Top = 40
        Width = 121
        Height = 21
        Hint = 'The depth of the potential hole, K'
        Alignment = taLeftJustify
        BevelInner = bvLowered
        BevelOuter = bvLowered
        Color = clWindow
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clGray
        Font.Height = -11
        Font.Name = 'MS Sans Serif'
        Font.Style = []
        ParentFont = False
        ParentShowHint = False
        ShowHint = True
        TabOrder = 2
      end
      object Potential_DEdit: TPanel
        Left = 130
        Top = 70
        Width = 121
        Height = 21
        Hint = 'The distance to the minimum of the potential'
        Alignment = taLeftJustify
        BevelInner = bvLowered
        BevelOuter = bvLowered
        Color = clWindow
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clGray
        Font.Height = -11
        Font.Name = 'MS Sans Serif'
        Font.Style = []
        ParentFont = False
        ParentShowHint = False
        ShowHint = True
        TabOrder = 3
      end
      object CellSize_2ndLabel: TPanel
        Left = 210
        Top = 100
        Width = 31
        Height = 21
        BevelOuter = bvNone
        Caption = 'by'
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -11
        Font.Name = 'MS Sans Serif'
        Font.Style = [fsBold]
        ParentFont = False
        TabOrder = 4
      end
      object CellsizeYEdit: TPanel
        Left = 250
        Top = 100
        Width = 71
        Height = 21
        Hint = 'Cell height'
        Alignment = taLeftJustify
        BevelInner = bvLowered
        BevelOuter = bvLowered
        Color = clWindow
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clGray
        Font.Height = -11
        Font.Name = 'MS Sans Serif'
        Font.Style = []
        ParentFont = False
        ParentShowHint = False
        ShowHint = True
        TabOrder = 5
      end
      object CellsizeXEdit: TPanel
        Left = 130
        Top = 100
        Width = 71
        Height = 21
        Hint = 'Cell width'
        Alignment = taLeftJustify
        BevelInner = bvLowered
        BevelOuter = bvLowered
        Color = clWindow
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clGray
        Font.Height = -11
        Font.Name = 'MS Sans Serif'
        Font.Style = []
        ParentFont = False
        ParentShowHint = False
        ShowHint = True
        TabOrder = 6
      end
      object CellSizeLabel: TPanel
        Left = 0
        Top = 100
        Width = 121
        Height = 21
        BevelOuter = bvNone
        Caption = 'Cell size'
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -11
        Font.Name = 'MS Sans Serif'
        Font.Style = [fsBold]
        ParentFont = False
        TabOrder = 7
      end
      object SPFLabel: TPanel
        Left = 0
        Top = 160
        Width = 121
        Height = 21
        BevelOuter = bvNone
        Caption = 'Steps per frame'
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -11
        Font.Name = 'MS Sans Serif'
        Font.Style = [fsBold]
        ParentFont = False
        TabOrder = 8
      end
      object SPFEdit: TEdit
        Left = 130
        Top = 160
        Width = 121
        Height = 21
        Hint = 'Steps per frame'
        AutoSelect = False
        ParentShowHint = False
        ShowHint = True
        TabOrder = 9
      end
      object SPFApplyButton: TButton
        Left = 260
        Top = 160
        Width = 61
        Height = 21
        Hint = 'Apply new steps per frame value'
        Caption = 'Apply'
        ParentShowHint = False
        ShowHint = True
        TabOrder = 10
        OnClick = SPFApplyButtonClick
      end
      object ApplyTimeStepButton: TButton
        Left = 260
        Top = 130
        Width = 61
        Height = 21
        Hint = 'Apply new simulation time step'
        Caption = 'Apply'
        ParentShowHint = False
        ShowHint = True
        TabOrder = 11
        OnClick = ApplyTimeStepButtonClick
      end
      object TimeStepEdit: TEdit
        Left = 130
        Top = 130
        Width = 121
        Height = 21
        Hint = 'Simulation time step'
        AutoSelect = False
        ParentShowHint = False
        ShowHint = True
        TabOrder = 12
      end
      object TimeStepLabel: TPanel
        Left = 0
        Top = 130
        Width = 121
        Height = 21
        BevelOuter = bvNone
        Caption = 'Simulation time step'
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -11
        Font.Name = 'MS Sans Serif'
        Font.Style = [fsBold]
        ParentFont = False
        TabOrder = 13
      end
      object BoundaryLabel: TPanel
        Left = 0
        Top = 220
        Width = 121
        Height = 31
        BevelOuter = bvNone
        Caption = 'Boundary'
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -11
        Font.Name = 'MS Sans Serif'
        Font.Style = [fsBold]
        ParentFont = False
        TabOrder = 14
      end
      object ParticlesLabel: TPanel
        Left = 0
        Top = 10
        Width = 121
        Height = 21
        BevelOuter = bvNone
        Caption = 'Number of particles'
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -11
        Font.Name = 'MS Sans Serif'
        Font.Style = [fsBold]
        ParentFont = False
        TabOrder = 15
      end
      object ParticlesView: TPanel
        Left = 130
        Top = 10
        Width = 121
        Height = 21
        Hint = 'Number of moleculs in the simulation'
        Alignment = taLeftJustify
        BevelInner = bvLowered
        BevelOuter = bvLowered
        Color = clWindow
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clGray
        Font.Height = -11
        Font.Name = 'MS Sans Serif'
        Font.Style = []
        ParentFont = False
        ParentShowHint = False
        ShowHint = True
        TabOrder = 16
      end
      object InitialKLabel: TPanel
        Left = 0
        Top = 190
        Width = 121
        Height = 21
        BevelOuter = bvNone
        Caption = 'Initial kinetic/atom, K'
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -11
        Font.Name = 'MS Sans Serif'
        Font.Style = [fsBold]
        ParentFont = False
        TabOrder = 17
      end
      object InitialKEdit: TEdit
        Left = 130
        Top = 190
        Width = 121
        Height = 21
        Hint = 'Initial kinetic energy per atom, K'
        AutoSelect = False
        ParentShowHint = False
        ShowHint = True
        TabOrder = 18
      end
      object InitialKApplyButton: TButton
        Left = 260
        Top = 190
        Width = 61
        Height = 21
        Hint = 'Apply new initial kinetic energy per atom value'
        Caption = 'Apply'
        ParentShowHint = False
        ShowHint = True
        TabOrder = 19
        OnClick = InitialKApplyButtonClick
      end
    end
    object StatisticsTab: TTabSheet
      Caption = 'Statistics'
      ImageIndex = 1
      object TotalTimeLabel: TPanel
        Left = 0
        Top = 10
        Width = 111
        Height = 21
        BevelOuter = bvNone
        Caption = 'Total time'
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -11
        Font.Name = 'MS Sans Serif'
        Font.Style = [fsBold]
        ParentFont = False
        TabOrder = 0
      end
      object TotalTimeView: TPanel
        Left = 120
        Top = 10
        Width = 201
        Height = 21
        Hint = 'Current total time in the simulation'
        Alignment = taLeftJustify
        BevelInner = bvLowered
        BevelOuter = bvLowered
        Color = clWindow
        Font.Charset = ANSI_CHARSET
        Font.Color = clGray
        Font.Height = -11
        Font.Name = 'Courier New'
        Font.Style = []
        ParentFont = False
        ParentShowHint = False
        ShowHint = True
        TabOrder = 1
      end
      object TotalEnergyLabel: TPanel
        Left = 0
        Top = 40
        Width = 111
        Height = 21
        BevelOuter = bvNone
        Caption = 'Total energy, K'
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -11
        Font.Name = 'MS Sans Serif'
        Font.Style = [fsBold]
        ParentFont = False
        TabOrder = 2
      end
      object TotalEnergyView: TPanel
        Left = 120
        Top = 40
        Width = 201
        Height = 21
        Hint = 'Mean energy per atom, K'
        Alignment = taLeftJustify
        BevelInner = bvLowered
        BevelOuter = bvLowered
        Color = clWindow
        Font.Charset = ANSI_CHARSET
        Font.Color = clGray
        Font.Height = -11
        Font.Name = 'Courier New'
        Font.Style = []
        ParentFont = False
        ParentShowHint = False
        ShowHint = True
        TabOrder = 3
      end
      object KEnergyLabel: TPanel
        Left = 0
        Top = 100
        Width = 111
        Height = 21
        BevelOuter = bvNone
        Caption = 'Kinetic energy, K'
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -11
        Font.Name = 'MS Sans Serif'
        Font.Style = [fsBold]
        ParentFont = False
        TabOrder = 4
      end
      object MinKEnergyView: TPanel
        Left = 120
        Top = 100
        Width = 61
        Height = 21
        Hint = 'Minimal kinetic energy of a particle, K'
        Alignment = taLeftJustify
        BevelInner = bvLowered
        BevelOuter = bvLowered
        Color = clWindow
        Font.Charset = ANSI_CHARSET
        Font.Color = clGray
        Font.Height = -11
        Font.Name = 'Courier New'
        Font.Style = []
        ParentFont = False
        ParentShowHint = False
        ShowHint = True
        TabOrder = 5
      end
      object MaxKEnergyView: TPanel
        Left = 260
        Top = 100
        Width = 61
        Height = 21
        Hint = 'Maximal kinetic energy of a particle, K'
        Alignment = taLeftJustify
        BevelInner = bvLowered
        BevelOuter = bvLowered
        Color = clWindow
        Font.Charset = ANSI_CHARSET
        Font.Color = clGray
        Font.Height = -11
        Font.Name = 'Courier New'
        Font.Style = []
        ParentFont = False
        ParentShowHint = False
        ShowHint = True
        TabOrder = 6
      end
      object MinUEnergyLabel: TPanel
        Left = 0
        Top = 130
        Width = 111
        Height = 21
        BevelOuter = bvNone
        Caption = 'Potential energy, K'
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -11
        Font.Name = 'MS Sans Serif'
        Font.Style = [fsBold]
        ParentFont = False
        TabOrder = 7
      end
      object MaxUEnergyView: TPanel
        Left = 260
        Top = 130
        Width = 61
        Height = 21
        Hint = 'Minimal potential energy of a particle, K'
        Alignment = taLeftJustify
        BevelInner = bvLowered
        BevelOuter = bvLowered
        Color = clWindow
        Font.Charset = ANSI_CHARSET
        Font.Color = clGray
        Font.Height = -11
        Font.Name = 'Courier New'
        Font.Style = []
        ParentFont = False
        ParentShowHint = False
        ShowHint = True
        TabOrder = 8
      end
      object MinUEnergyView: TPanel
        Left = 120
        Top = 130
        Width = 61
        Height = 21
        Hint = 'Maximal potential energy of a particle, K'
        Alignment = taLeftJustify
        BevelInner = bvLowered
        BevelOuter = bvLowered
        Color = clWindow
        Font.Charset = ANSI_CHARSET
        Font.Color = clGray
        Font.Height = -11
        Font.Name = 'Courier New'
        Font.Style = []
        ParentFont = False
        ParentShowHint = False
        ShowHint = True
        TabOrder = 9
      end
      object ErrorLabel: TPanel
        Left = 0
        Top = 70
        Width = 111
        Height = 21
        BevelOuter = bvNone
        Caption = 'Energy error'
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -11
        Font.Name = 'MS Sans Serif'
        Font.Style = [fsBold]
        ParentFont = False
        TabOrder = 10
      end
      object ErrorView: TPanel
        Left = 120
        Top = 70
        Width = 201
        Height = 21
        Hint = '(energy - initial_energy) / initial_energy'
        Alignment = taLeftJustify
        BevelInner = bvLowered
        BevelOuter = bvLowered
        Color = clWindow
        Font.Charset = ANSI_CHARSET
        Font.Color = clGray
        Font.Height = -11
        Font.Name = 'Courier New'
        Font.Style = []
        ParentFont = False
        ParentShowHint = False
        ShowHint = True
        TabOrder = 11
      end
      object DLabel: TPanel
        Left = 0
        Top = 190
        Width = 111
        Height = 21
        BevelOuter = bvNone
        Caption = 'Diffusion index'
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -11
        Font.Name = 'MS Sans Serif'
        Font.Style = [fsBold]
        ParentFont = False
        TabOrder = 12
      end
      object DView: TPanel
        Left = 120
        Top = 190
        Width = 201
        Height = 21
        Hint = '<x^2>/2t'
        Alignment = taLeftJustify
        BevelInner = bvLowered
        BevelOuter = bvLowered
        Color = clWindow
        Font.Charset = ANSI_CHARSET
        Font.Color = clGray
        Font.Height = -11
        Font.Name = 'Courier New'
        Font.Style = []
        ParentFont = False
        ParentShowHint = False
        ShowHint = True
        TabOrder = 13
      end
      object PressureLabel: TPanel
        Left = 0
        Top = 220
        Width = 111
        Height = 21
        BevelOuter = bvNone
        Caption = 'Average pressure'
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -11
        Font.Name = 'MS Sans Serif'
        Font.Style = [fsBold]
        ParentFont = False
        TabOrder = 14
      end
      object XPressureView: TPanel
        Left = 120
        Top = 220
        Width = 91
        Height = 21
        Hint = 'x pressure'
        Alignment = taLeftJustify
        BevelInner = bvLowered
        BevelOuter = bvLowered
        Color = clWindow
        Font.Charset = ANSI_CHARSET
        Font.Color = clGray
        Font.Height = -11
        Font.Name = 'Courier New'
        Font.Style = []
        ParentFont = False
        ParentShowHint = False
        ShowHint = True
        TabOrder = 15
      end
      object YPressureView: TPanel
        Left = 230
        Top = 220
        Width = 91
        Height = 21
        Hint = 'y pressure'
        Alignment = taLeftJustify
        BevelInner = bvLowered
        BevelOuter = bvLowered
        Color = clWindow
        Font.Charset = ANSI_CHARSET
        Font.Color = clGray
        Font.Height = -11
        Font.Name = 'Courier New'
        Font.Style = []
        ParentFont = False
        ParentShowHint = False
        ShowHint = True
        TabOrder = 16
      end
      object DR2Label: TPanel
        Left = 0
        Top = 160
        Width = 111
        Height = 21
        BevelOuter = bvNone
        Caption = 'Mean (R-R0)^2'
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -11
        Font.Name = 'MS Sans Serif'
        Font.Style = [fsBold]
        ParentFont = False
        TabOrder = 17
      end
      object DR2View: TPanel
        Left = 120
        Top = 160
        Width = 201
        Height = 21
        Hint = '<x^2>'
        Alignment = taLeftJustify
        BevelInner = bvLowered
        BevelOuter = bvLowered
        Color = clWindow
        Font.Charset = ANSI_CHARSET
        Font.Color = clGray
        Font.Height = -11
        Font.Name = 'Courier New'
        Font.Style = []
        ParentFont = False
        ParentShowHint = False
        ShowHint = True
        TabOrder = 18
      end
      object AveKEnergyView: TPanel
        Left = 190
        Top = 100
        Width = 61
        Height = 21
        Hint = 'Average kinetic energy of a particle, K'
        Alignment = taLeftJustify
        BevelInner = bvLowered
        BevelOuter = bvLowered
        Color = clWindow
        Font.Charset = ANSI_CHARSET
        Font.Color = clGray
        Font.Height = -11
        Font.Name = 'Courier New'
        Font.Style = []
        ParentFont = False
        ParentShowHint = False
        ShowHint = True
        TabOrder = 19
      end
      object AveUEnergyView: TPanel
        Left = 190
        Top = 130
        Width = 61
        Height = 21
        Hint = 'Average potential energy of a particle, K'
        Alignment = taLeftJustify
        BevelInner = bvLowered
        BevelOuter = bvLowered
        Color = clWindow
        Font.Charset = ANSI_CHARSET
        Font.Color = clGray
        Font.Height = -11
        Font.Name = 'Courier New'
        Font.Style = []
        ParentFont = False
        ParentShowHint = False
        ShowHint = True
        TabOrder = 20
      end
    end
    object GRTab: TTabSheet
      Caption = 'Radial distribution'
      ImageIndex = 2
      object GRPanel: TPanel
        Left = 0
        Top = 0
        Width = 331
        Height = 251
        Hint = 'g(r), the radial distribution'
        BevelOuter = bvLowered
        Color = clWhite
        ParentShowHint = False
        ShowHint = True
        TabOrder = 0
        object GRView: TImage
          Left = 0
          Top = 0
          Width = 331
          Height = 251
          HelpContext = 1
          ParentShowHint = False
          ShowHint = False
        end
      end
    end
    object KTab: TTabSheet
      Caption = 'Kinetic energy'
      ImageIndex = 3
      object KEnergyPanel: TPanel
        Left = 0
        Top = 0
        Width = 331
        Height = 251
        Hint = 'Kinetic energy distribution across particles'
        BevelOuter = bvLowered
        Color = clWhite
        ParentShowHint = False
        ShowHint = True
        TabOrder = 0
        object KEnergyView: TImage
          Left = 0
          Top = 0
          Width = 331
          Height = 251
          HelpContext = 1
          ParentShowHint = False
          ShowHint = False
        end
      end
    end
    object UTab: TTabSheet
      Caption = 'Potential energy'
      ImageIndex = 4
      object UEnergyPanel: TPanel
        Left = 0
        Top = 0
        Width = 331
        Height = 251
        Hint = 'Potential energy distribution across perticles'
        BevelOuter = bvLowered
        Color = clWhite
        ParentShowHint = False
        ShowHint = True
        TabOrder = 0
        object UEnergyView: TImage
          Left = 0
          Top = 0
          Width = 331
          Height = 251
          HelpContext = 1
          ParentShowHint = False
          ShowHint = False
        end
      end
    end
    object ETab: TTabSheet
      Caption = 'Full energy'
      ImageIndex = 5
      object EnergyPanel: TPanel
        Left = 0
        Top = 0
        Width = 331
        Height = 251
        Hint = 'Full energy distribution across perticles'
        BevelOuter = bvLowered
        Color = clWhite
        ParentShowHint = False
        ShowHint = True
        TabOrder = 0
        object EnergyView: TImage
          Left = 0
          Top = 0
          Width = 331
          Height = 251
          HelpContext = 1
          ParentShowHint = False
          ShowHint = False
        end
      end
    end
  end
end
