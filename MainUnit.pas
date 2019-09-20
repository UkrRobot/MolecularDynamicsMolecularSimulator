unit MainUnit;

interface

uses
  Windows, Messages, SysUtils, Variants, Classes, Graphics, Controls, Forms,
  Dialogs, StdCtrls, ExtCtrls, Buttons, Math, ComCtrls;

const
  N = 100; { total number of atoms }
  cellsize_x = 2;  //cellsize true= dtrue*10
  cellsize_y = 2;
  { potential is a/r^12-b/r^6 = e((d/r)^12-2(d/r)^6) }
  { it's minimum is at (2*a/b)^(1/6) = d }
  { the minimal value is -b^2/(4*a) = -e }
  potential_d = 0.2; // potential_a = 100*d^6 если  dtrue=3.4e-10 то dt=dt*??????
  potential_e = 100;   //d==sigma      potenciall==force
  potential_a = potential_e*potential_d*potential_d*potential_d*potential_d*
                potential_d*potential_d*potential_d*potential_d*potential_d*
                potential_d*potential_d*potential_d;
  potential_b = 2*potential_e*potential_d*potential_d*potential_d*potential_d*
                potential_d*potential_d;
  border_k = 2*potential_e / 0.004;

  initial_x1 = 0.0; initial_x2 = 2.0;
  initial_y1 = 0.0; initial_y2 = 2.0;
  initial_dt = 0.000005;
  initial_kinetic = potential_e;


  Graph_N = 300;

  { for g(r) calculations }
  GR_max_R = potential_d * 8;
  GR_N = Graph_N;
  GR_step = GR_max_R/(GR_N-1);
  GR_shift = 0;
  GR_strip_width = potential_d/4.0;
  { for n(e) calculations }
  E_points = Graph_N;
  E_slices = 10;

const
  { potential modes }
  Mode_Fixed = 1;
  Mode_Periodic = 2;

type
  TMainForm = class(TForm)
    ImagePanel: TPanel;
    AtomView: TImage;
    StartButton: TButton;
    ResetButton: TButton;
    QuitButton: TButton;
    OneFrameButton: TButton;
    AdditionalPages: TPageControl;
    ParametersTab: TTabSheet;
    StatisticsTab: TTabSheet;
    GRTab: TTabSheet;
    Potential_ELabel: TPanel;
    Potential_DLabel: TPanel;
    Potential_EEdit: TPanel;
    Potential_DEdit: TPanel;
    CellSize_2ndLabel: TPanel;
    CellsizeYEdit: TPanel;
    CellsizeXEdit: TPanel;
    CellSizeLabel: TPanel;
    SPFLabel: TPanel;
    SPFEdit: TEdit;
    SPFApplyButton: TButton;
    ApplyTimeStepButton: TButton;
    TimeStepEdit: TEdit;
    TimeStepLabel: TPanel;
    PeriodicBorderSwitch: TSpeedButton;
    FixedBorderSwitch: TSpeedButton;
    BoundaryLabel: TPanel;
    GRPanel: TPanel;
    GRView: TImage;
    TotalTimeLabel: TPanel;
    TotalTimeView: TPanel;
    TotalEnergyLabel: TPanel;
    TotalEnergyView: TPanel;
    KEnergyLabel: TPanel;
    MinKEnergyView: TPanel;
    MaxKEnergyView: TPanel;
    KTab: TTabSheet;
    KEnergyPanel: TPanel;
    KEnergyView: TImage;
    ParticlesLabel: TPanel;
    ParticlesView: TPanel;
    MinUEnergyLabel: TPanel;
    MaxUEnergyView: TPanel;
    MinUEnergyView: TPanel;
    UTab: TTabSheet;
    UEnergyPanel: TPanel;
    UEnergyView: TImage;
    ETab: TTabSheet;
    EnergyPanel: TPanel;
    EnergyView: TImage;
    ErrorLabel: TPanel;
    ErrorView: TPanel;
    DLabel: TPanel;
    DView: TPanel;
    InitialKLabel: TPanel;
    InitialKEdit: TEdit;
    InitialKApplyButton: TButton;
    PressureLabel: TPanel;
    XPressureView: TPanel;
    YPressureView: TPanel;
    DR2Label: TPanel;
    DR2View: TPanel;
    AveKEnergyView: TPanel;
    AveUEnergyView: TPanel;
    procedure StartButtonClick(Sender: TObject);
    procedure FormCreate(Sender: TObject);
    procedure FormDestroy(Sender: TObject);
    procedure ResetButtonClick(Sender: TObject);
    procedure QuitButtonClick(Sender: TObject);
    procedure FixedBorderSwitchClick(Sender: TObject);
    procedure PeriodicBorderSwitchClick(Sender: TObject);
    procedure SPFApplyButtonClick(Sender: TObject);
    procedure OneFrameButtonClick(Sender: TObject);
    procedure AdditionalPagesChange(Sender: TObject);
    procedure ApplyTimeStepButtonClick(Sender: TObject);
    procedure InitialKApplyButtonClick(Sender: TObject);
  private
  { Private declarations }
    Running, StopFlag: boolean;
    Runner: TThread;
    time_step, initial_K: double;
    SPF: integer;
    procedure StopRunning;
    procedure StartRunning;
    procedure RedrawImage;
    procedure DrawAtoms(Canvas: TCanvas);
    procedure UpdateStatistics;
  public
    { Public declarations }
  end;

  TRunner = class(TThread)
  public
    Form: TMainForm;
  protected
    procedure Execute; override;
  end;

  Atom = record
    x, y: double;
    vx, vy: double;
  end;
  Energy = record
    u, k: double;
  end;

  array_of_double = array of double;

  Atom_Array = array [1..N] of Atom;
  Energy_Array = array [1..N] of Energy;
  Graph_Array = array [0..Graph_N-1] of double;

var
  MainForm: TMainForm;
  initial_atoms, atoms: Atom_Array;
  initial_energy, impulse_x, impulse_y, total_time: double;

  k1, k2, k3, k4, k5, k6, temp_atoms: Atom_Array;
  EA: Energy_Array;
  GR: Graph_Array;
  mode: integer;

implementation

{$R *.dfm}

function rand_double(min, max: double): double;
const granularity = 30000;
begin
  rand_double := min+((max-min)*random(granularity))/(granularity-1);
end;

{ ----------------------- Basic Force & Potential ----------------------- }

{ a/r^12-b/r^6 = (a/r^6-b)/r^6 }
function potential_direct(dx, dy: double): double;
var d2,d6: double;
begin
  d2 := 1.0/(sqr(dx)+sqr(dy));
  d6 := d2*d2*d2;
  potential_direct := (potential_a*d6 - potential_b)*d6;
end;

function border_potential(x, y: double): double;
var p: double;
begin
  p := 0;
  if x > cellsize_x then p := p + border_k*(x-cellsize_x)*(x-cellsize_x)/2.0
  else if x < 0 then p := p + border_k*x*x/2.0;
  if y > cellsize_y then p := p + border_k*(y-cellsize_y)*(y-cellsize_y)/2.0
  else if y < 0 then p := p + border_k*y*y/2.0;
  border_potential := p;
end;

{ 12*a/r^13-6*b/r^7 = (12*a/r^6-6*b)/r^7 }
procedure force_direct(dx, dy: double; var fx, fy: double);
var d2,d6,f: double;
begin
  d2 := 1.0/(sqr(dx)+sqr(dy));
  d6 := d2*d2*d2;
  f := (12*potential_a*d6 - 6*potential_b)*d6*d2;  //force
  fx := f*dx;
  fy := f*dy;
end;

procedure border_force(x, y: double; var fx, fy: double);
begin
  if x > cellsize_x then fx := -border_k*(x-cellsize_x)
  else if x < 0 then fx := -border_k*x
  else fx := 0;
  if y > cellsize_y then fy := -border_k*(y-cellsize_y)
  else if y < 0 then fy := -border_k*y
  else fy := 0;
end;

{ the cutoff distance is concidered to be less than half of cell size }
function potential_periodic(dx, dy: double): double;
begin
   dx := abs(dx);
   dy := abs(dy);
   if dx > cellsize_x*0.5 then dx := dx - cellsize_x;
   if dy > cellsize_y*0.5 then dy := dy - cellsize_y;
   potential_periodic := potential_direct(dx, dy);
end;

{ the cutoff distance is concidered to be less than half of cell size }
procedure force_periodic(dx, dy: double; var fx, fy: double);
var ddx, ddy: double;
begin
   ddx := abs(dx);
   ddy := abs(dy);
   if ddx > cellsize_x*0.5 then ddx := ddx - cellsize_x;
   if ddy > cellsize_y*0.5 then ddy := ddy - cellsize_y;

   force_direct(ddx, ddy, fx, fy);
   if dx < 0 then fx := -fx;
   if dy < 0 then fy := -fy;
end;

{ ----------------------- Forces & Potentials ----------------------- }

procedure calc_derivatives(a: Atom_Array; var f: Atom_Array);
var i, j: integer;
    fx, fy: double;      // derivatives=приращения
begin
  for i := 1 to N do
  begin
    f[i].x := a[i].vx;   //вспомогательній массив приращения координат за 1
    f[i].y := a[i].vy;
    f[i].vx := 0;        //приращения скоросте за 1 времени  1=????
    f[i].vy := 0;
  end;
  if mode = Mode_Fixed then
  begin { fixed border }
    for i := 1 to N-1 do
    for j := i+1 to N do
    begin
      force_direct(a[i].x-a[j].x, a[i].y-a[j].y, fx, fy);
      f[i].vx := f[i].vx + fx;    //m=1 ???? where dt???    fx==ax
      f[i].vy := f[i].vy + fy;      //vx=ax????
      f[j].vx := f[j].vx - fx;
      f[j].vy := f[j].vy - fy;
    end;
    for i := 1 to N do
    begin
      border_force(a[i].x, a[i].y, fx, fy);
      f[i].vx := f[i].vx + fx;
      f[i].vy := f[i].vy + fy;
    end;
  end else begin  { periodic border }
    for i := 1 to N-1 do
    for j := i+1 to N do
    begin
      force_periodic(a[i].x-a[j].x, a[i].y-a[j].y, fx, fy);
      f[i].vx := f[i].vx + fx;
      f[i].vy := f[i].vy + fy;
      f[j].vx := f[j].vx - fx;
      f[j].vy := f[j].vy - fy;
    end;
  end;
end;

procedure calc_energy(a: Atom_Array; var EA: Energy_Array;
                      var kmin, kmax, kave: double;
                      var umin, umax, uave: double;
                      var emin, emax, eave, etotal: double);
var
  i, j: integer;
  e: double;
begin
  etotal := 0;

  for i := 1 to N do
  begin
    EA[i].k := (sqr(a[i].vx)+sqr(a[i].vy))*0.5;
    EA[i].u := 0;
    etotal := etotal + EA[i].k;
  end;

  if mode = Mode_Fixed then
  begin { fixed border }
    for i := 1 to N-1 do
    for j := i+1 to N do
    begin
      e := potential_direct(a[i].x-a[j].x, a[i].y-a[j].y);
      EA[i].u := EA[i].u + e;
      EA[j].u := EA[j].u + e;
      etotal := etotal + e;
    end;
    for i := 1 to N do
    begin
      e := border_potential(a[i].x, a[i].y);
      EA[i].u := EA[i].u + e;
      etotal := etotal + e;
    end;
  end else begin  { periodic border }
    for i := 1 to N-1 do
    for j := i+1 to N do
    begin
      e := potential_periodic(a[i].x-a[j].x, a[i].y-a[j].y);
      EA[i].u := EA[i].u + e;
      EA[j].u := EA[j].u + e;
      etotal := etotal + e;
    end;
  end;

  kmin := +1e20; kmax := -1e20; kave := 0;
  umin := +1e20; umax := -1e20; uave := 0;
  emin := +1e20; emax := -1e20; eave := 0;

  for i := 1 to N do
  begin
    e := EA[i].k;
    if e > kmax then kmax := e
    else if e < kmin then kmin := e;
    kave := kave + e;

    e := EA[i].u;
    if e > umax then umax := e
    else if e < umin then umin := e;
    uave := uave + e;

    e := EA[i].u + EA[i].k;
    if e > emax then emax := e
    else if e < emin then emin := e;
    eave := eave + e;
  end;
  kave := kave / N;
  uave := uave / N;
  eave := eave / N;
end;

{ ----------------------- Model Calculations ----------------------- }

procedure rk5step(a: Atom_Array; var out_atoms: Atom_Array; dt: double);
var i: integer;
    fx, fy: double;
const  b21 = 1.0/5.0;     //5pixels????
  b31 = 3.0/40.0;
  b32 = 9.0/40.0;
  b41 = 3.0/10.0;
  b42 = -9.0/10.0;
  b43 = 6.0/5.0;
  b51 = -11.0/54.0;
  b52 = 5.0/2.0;
  b53 = -70.0/27.0;
  b54 = 35.0/27.0;
  b61 = 1631.0/55296.0;
  b62 = 175.0/512.0;
  b63 = 575.0/13824.0;
  b64 = 44275.0/110592.0;
  b65 = 253.0/4096.0;
  c1 = 37.0/378.0;
  c2 = 0.0;
  c3 = 250.0/621.0;
  c4 = 125.0/594.0;
  c5 = 0.0;
  c6 = 512.0/1771.0;
begin
  calc_derivatives(a, k1);    //приращений
  for i := 1 to N do
  begin
    temp_atoms[i].x := a[i].x + (b21*k1[i].x)*dt;     //
    temp_atoms[i].y := a[i].y + (b21*k1[i].y)*dt;
    temp_atoms[i].vx := a[i].vx + (b21*k1[i].vx)*dt;
    temp_atoms[i].vy := a[i].vy + (b21*k1[i].vy)*dt;
  end;
  calc_derivatives(temp_atoms, k2);
  for i := 1 to N do
  begin
    temp_atoms[i].x := a[i].x + (b31*k1[i].x + b32*k2[i].x)*dt;
    temp_atoms[i].y := a[i].y + (b31*k1[i].y + b32*k2[i].y)*dt;
    temp_atoms[i].vx := a[i].vx + (b31*k1[i].vx + b32*k2[i].vx)*dt;
    temp_atoms[i].vy := a[i].vy + (b31*k1[i].vy + b32*k2[i].vy)*dt;
  end;
  calc_derivatives(temp_atoms, k3);
  for i := 1 to N do
  begin
    temp_atoms[i].x := a[i].x + (b41*k1[i].x + b42*k2[i].x + b43*k3[i].x)*dt;
    temp_atoms[i].y := a[i].y + (b41*k1[i].y + b42*k2[i].y + b43*k3[i].y)*dt;
    temp_atoms[i].vx := a[i].vx + (b41*k1[i].vx + b42*k2[i].vx + b43*k3[i].vx)*dt;
    temp_atoms[i].vy := a[i].vy + (b41*k1[i].vy + b42*k2[i].vy + b43*k3[i].vy)*dt;
  end;
  calc_derivatives(temp_atoms, k4);
  for i := 1 to N do
  begin
    temp_atoms[i].x := a[i].x + (b51*k1[i].x + b52*k2[i].x + b53*k3[i].x + b54*k4[i].x)*dt;
    temp_atoms[i].y := a[i].y + (b51*k1[i].y + b52*k2[i].y + b53*k3[i].y + b54*k4[i].y)*dt;
    temp_atoms[i].vx := a[i].vx + (b51*k1[i].vx + b52*k2[i].vx + b53*k3[i].vx + b54*k4[i].vx)*dt;
    temp_atoms[i].vy := a[i].vy + (b51*k1[i].vy + b52*k2[i].vy + b53*k3[i].vy + b54*k4[i].vy)*dt;
  end;
  calc_derivatives(temp_atoms, k5);
  for i := 1 to N do
  begin
    temp_atoms[i].x := a[i].x + (b61*k1[i].x + b62*k2[i].x + b63*k3[i].x + b64*k4[i].x + b65*k5[i].x)*dt;
    temp_atoms[i].y := a[i].y + (b61*k1[i].y + b62*k2[i].y + b63*k3[i].y + b64*k4[i].y + b65*k5[i].y)*dt;
    temp_atoms[i].vx := a[i].vx + (b61*k1[i].vx + b62*k2[i].vx + b63*k3[i].vx + b64*k4[i].vx + b65*k5[i].vx)*dt;
    temp_atoms[i].vy := a[i].vy + (b61*k1[i].vy + b62*k2[i].vy + b63*k3[i].vy + b64*k4[i].vy + b65*k5[i].vy)*dt;
  end;
  calc_derivatives(temp_atoms, k6);
  if mode = Mode_Fixed then         //упругие фиксированные стенки
  begin
    for i := 1 to N do
    begin
      out_atoms[i].x := a[i].x + (c1*k1[i].x + c2*k2[i].x + c3*k3[i].x + c4*k4[i].x + c5*k5[i].x + c6*k6[i].x)*dt;
      out_atoms[i].y := a[i].y + (c1*k1[i].y + c2*k2[i].y + c3*k3[i].y + c4*k4[i].y + c5*k5[i].y + c6*k6[i].y)*dt;
      out_atoms[i].vx := a[i].vx + (c1*k1[i].vx + c2*k2[i].vx + c3*k3[i].vx + c4*k4[i].vx + c5*k5[i].vx + c6*k6[i].vx)*dt;
      out_atoms[i].vy := a[i].vy + (c1*k1[i].vy + c2*k2[i].vy + c3*k3[i].vy + c4*k4[i].vy + c5*k5[i].vy + c6*k6[i].vy)*dt;

      border_force(out_atoms[i].x, out_atoms[i].y, fx, fy);
      impulse_x := impulse_x + abs(fx);        // добавка к давлению на кожному кроци при упругом
      impulse_y := impulse_y + abs(fy);
    end;
  end else begin		// переодические граничные
    for i := 1 to N do
    begin
      out_atoms[i].x := a[i].x + (c1*k1[i].x + c2*k2[i].x + c3*k3[i].x + c4*k4[i].x + c5*k5[i].x + c6*k6[i].x)*dt;
      out_atoms[i].y := a[i].y + (c1*k1[i].y + c2*k2[i].y + c3*k3[i].y + c4*k4[i].y + c5*k5[i].y + c6*k6[i].y)*dt;
      out_atoms[i].vx := a[i].vx + (c1*k1[i].vx + c2*k2[i].vx + c3*k3[i].vx + c4*k4[i].vx + c5*k5[i].vx + c6*k6[i].vx)*dt;
      out_atoms[i].vy := a[i].vy + (c1*k1[i].vy + c2*k2[i].vy + c3*k3[i].vy + c4*k4[i].vy + c5*k5[i].vy + c6*k6[i].vy)*dt;

      if out_atoms[i].x > cellsize_x then
      begin
        out_atoms[i].x := out_atoms[i].x - cellsize_x;
        initial_atoms[i].x := initial_atoms[i].x - cellsize_x;
        impulse_x := impulse_x + abs(2*out_atoms[i].vx);
      end;
      if out_atoms[i].y > cellsize_y then
      begin
        out_atoms[i].y := out_atoms[i].y - cellsize_y;
        initial_atoms[i].y := initial_atoms[i].y - cellsize_y;
        impulse_y := impulse_y + abs(2*out_atoms[i].vy);
      end;
      if out_atoms[i].x < 0 then
      begin
        out_atoms[i].x := out_atoms[i].x + cellsize_x;
        initial_atoms[i].x := initial_atoms[i].x + cellsize_x;
        impulse_x := impulse_x + abs(2*out_atoms[i].vx);   // добавка к давлению после периодических
      end;
      if out_atoms[i].y < 0 then
      begin
        out_atoms[i].y := out_atoms[i].y + cellsize_y;
        initial_atoms[i].y := initial_atoms[i].y + cellsize_y;
        impulse_y := impulse_y + abs(2*out_atoms[i].vy);
      end;
    end;
  end;
end;

{ ----------------------- Model Interface ----------------------- }

procedure calc_GR(var GR: Graph_Array;
                  var gmin, gmax, xmin, xmax: double);
var
  i, j, k, kmin, kmax: integer;
  c: array [0.. Graph_N-1] of integer;
  rmin, rmax, rmin2, rmax2: double;
  s, d, d2, koef: double;
begin
  for k := 0 to Graph_N-1 do c[k] := 0;

  for i := 1 to N-1 do
  for j := i+1 to N do
  begin
    d2 := sqr(atoms[i].x-atoms[j].x) + sqr(atoms[i].y-atoms[j].y);
    d := sqrt(d2);
    kmin := floor((d - GR_shift - GR_strip_width*0.5)/GR_step)-1;
    kmax := ceil((d - GR_shift + GR_strip_width*0.5)/GR_step)+1;
    if kmin < 0 then kmin := 0;
    if kmax > GR_N-1 then kmax := GR_N-1;
    for k := kmin to kmax do
    begin
      rmin := GR_shift + k*GR_step - GR_strip_width*0.5;
      rmax := rmin + GR_strip_width;
      if rmin < 0 then rmin := 0;
      rmin2 := sqr(rmin);
      rmax2 := sqr(rmax);
      if (d2 <= rmax2) and (d2 >= rmin2) then c[k] := c[k] + 2;
    end;
  end;

  gmin := +1e20;
  gmax := -1e20;
  xmin := 0;
  xmax := GR_max_R/potential_d;

  koef := cellsize_x*cellsize_y/N/N;

  for k := 0 to Graph_N-1 do
  begin
    rmin := GR_shift + k*GR_step - GR_strip_width*0.5;
    rmax := rmin + GR_strip_width;
    if rmin < 0 then rmin := 0;
    rmin2 := sqr(rmin);
    rmax2 := sqr(rmax);
    s := Pi*(rmax2-rmin2);
    GR[k] := koef*c[k]/(s);
    if GR[k] < gmin then gmin := GR[k];
    if GR[k] > gmax then gmax := GR[k];
  end;
end;

procedure calc_nK(var nE: Graph_Array;
                  var nmin, nmax, e_min, e_max: double);
var
  EA: Energy_Array;
  kmin, kmax, kave, umin, umax, uave, emin, emax, eave, etotal: double;
  min_k, max_k, i, k: integer;
begin
  calc_energy(atoms, EA, kmin, kmax, kave, umin, umax, uave, emin, emax, eave, etotal);
  e_min := kmin;
  e_max := kmax;

  for i := 0 to Graph_N-1 do nE[i] := 0;

  nmin := 0;
  nmax := 0;
  for i := 1 to N do
  begin
    min_k := floor(Graph_N*((EA[i].k-e_min)/(e_max-e_min) - 0.5/E_Slices));
    max_k := ceil(Graph_N*((EA[i].k-e_min)/(e_max-e_min) + 0.5/E_Slices));
    if min_k < 0 then min_k := 0;
    if max_k < 0 then max_k := 0;
    if min_k > Graph_N-1 then min_k := Graph_N-1;
    if max_k > Graph_N-1 then max_k := Graph_N-1;
    for k := min_k to max_k do
    begin
      nE[k] := nE[k] + 1;
      if nE[k] > nmax then nmax := nE[k];
    end;
  end;
  for i := 0 to Graph_N-1 do nE[i] := 100*nE[i]/N;
  nmin := 100*nmin/N;
  nmax := 100*nmax/N;
end;

procedure calc_nU(var nE: Graph_Array;
                  var nmin, nmax, e_min, e_max: double);
var
  EA: Energy_Array;
  kmin, kmax, kave, umin, umax, uave, emin, emax, eave, etotal: double;
  min_k, max_k, i, k: integer;
begin
  calc_energy(atoms, EA, kmin, kmax, kave, umin, umax, uave, emin, emax, eave, etotal);
  e_min := umin;
  e_max := umax;

  for i := 0 to Graph_N-1 do nE[i] := 0;

  nmin := 0;
  nmax := 0;
  for i := 1 to N do
  begin
    min_k := floor(Graph_N*((EA[i].u-e_min)/(e_max-e_min) - 0.5/E_Slices));
    max_k := ceil(Graph_N*((EA[i].u-e_min)/(e_max-e_min) + 0.5/E_Slices));
    if min_k < 0 then min_k := 0;
    if max_k < 0 then max_k := 0;
    if min_k > Graph_N-1 then min_k := Graph_N-1;
    if max_k > Graph_N-1 then max_k := Graph_N-1;
    for k := min_k to max_k do
    begin
      nE[k] := nE[k] + 1;
      if nE[k] > nmax then nmax := nE[k];
    end;
  end;
  for i := 0 to Graph_N-1 do nE[i] := 100*nE[i]/N;
  nmin := 100*nmin/N;
  nmax := 100*nmax/N;
end;

procedure calc_nE(var nE: Graph_Array;
                  var nmin, nmax, e_min, e_max: double);
var
  EA: Energy_Array;
  kmin, kmax, kave, umin, umax, uave, emin, emax, eave, etotal: double;
  min_k, max_k, i, k: integer;
begin
  calc_energy(atoms, EA, kmin, kmax, kave, umin, umax, uave, emin, emax, eave, etotal);
  e_min := emin;
  e_max := emax;

  for i := 0 to Graph_N-1 do nE[i] := 0;

  nmin := 0;
  nmax := 0;
  for i := 1 to N do
  begin
    min_k := floor(Graph_N*((EA[i].u+EA[i].k-e_min)/(e_max-e_min) - 0.5/E_Slices));
    max_k := ceil(Graph_N*((EA[i].u+EA[i].k-e_min)/(e_max-e_min) + 0.5/E_Slices));
    if min_k < 0 then min_k := 0;
    if max_k < 0 then max_k := 0;
    if min_k > Graph_N-1 then min_k := Graph_N-1;
    if max_k > Graph_N-1 then max_k := Graph_N-1;
    for k := min_k to max_k do
    begin
      nE[k] := nE[k] + 1;
      if nE[k] > nmax then nmax := nE[k];
    end;
  end;
  for i := 0 to Graph_N-1 do nE[i] := 100*nE[i]/N;
  nmin := 100*nmin/N;
  nmax := 100*nmax/N;
end;

procedure make_step(steps: integer; dt: double);
var i: integer;
begin
  for i := 1 to steps do
    rk5step(atoms, atoms, dt);
  total_time := total_time + dt * steps;
end;

procedure calc_deltas(var ave_dr2, ave_vx2, ave_vy2: double);
var dr2, vx2, vy2: double;
    i: integer;
begin
  dr2 := 0;
  vx2 := 0;
  vy2 := 0;
  for i := 1 to N do
  begin
    dr2 := dr2 + sqr(atoms[i].x-initial_atoms[i].x)+ sqr(atoms[i].y-initial_atoms[i].y);
    vx2 := vx2 + atoms[i].vx*initial_atoms[i].vx;
    vy2 := vy2 + atoms[i].vy*initial_atoms[i].vy;
  end;
  ave_dr2 := dr2/N;
  ave_vx2 := vx2/N;
  ave_vy2 := vy2/N;
end;

procedure reset_model(x1, x2, y1, y2, KE: double);
var i, j, k, kx, ky: integer;
    dx, dy, sx, sy: double;
    v, phi, dv: double;
    kmin, kmax, kave, umin, umax, uave, emin, emax, eave: double;
begin
  kx := round(sqrt(N));
  ky := ceil(N/kx);
  dx := (x2-x1)/(kx);
  dy := (y2-y1)/(ky);
  sx := x1 + dx/2;
  sy := y1 + dy/2;
  dv := sqrt(2*KE);
  i := 0;
  j := 0;
  for k := 1 to N do
  begin
    atoms[k].x := sx + dx*i;
    atoms[k].y := sy + dy*j;
    v := (rand_double(dv*0.9, dv*1.1)+rand_double(dv*0.9, dv*1.1)+
          rand_double(dv*0.9, dv*1.1)+rand_double(dv*0.9, dv*1.1))/4.0;
    phi := rand_double(0, 2*Pi);
    atoms[k].vx := v*cos(phi);
    atoms[k].vy := v*sin(phi);

    initial_atoms[k].x := atoms[k].x;
    initial_atoms[k].y := atoms[k].y;
    initial_atoms[k].vx := atoms[k].vx;
    initial_atoms[k].vy := atoms[k].vy;

    j := j + 1;
    if j >= ky then
    begin
      j := 0;
      i := i + 1;
    end;
  end;
  calc_energy(atoms, EA, kmin, kmax, kave, umin, umax, uave, emin, emax, eave, initial_energy);
  impulse_x := 0;
  impulse_y := 0;
  total_time := 0;
end;

{ ----------------------- Graph Drawing ----------------------- }

procedure guess_grid(min, max: double; var step, vmin, vmax: double);
var base: double;
begin
  base := power(10, floor(log10(abs(max-min))));
  vmin := floor(min/base)*base;
  vmax := ceil(max/base)*base;
  if abs(vmax-vmin) <= 1.2 * base then
  begin
    base := base / 5;
    vmin := floor(min/base)*base;
    vmax := ceil(max/base)*base;
  end;
  if abs(vmax-vmin) <= 3.0 * base then
  begin
    base := base / 2;
    vmin := floor(min/base)*base;
    vmax := ceil(max/base)*base;
  end;
  step := base;
end;

procedure DrawGrid(Canvas: TCanvas;
                    xtitle, ytitle: string;
                    var xmin, ymin, xmax, ymax: integer;
                    xvalmin, yvalmin, xvalmax, yvalmax: double);
var
  v, step: double;
  xvmin, xvmax, yvmin, yvmax: double;
  c: integer;
  s: string;
begin
  with Canvas do
  begin
    FillRect(ClipRect);
    MoveTo(xmin, ymax); LineTo(xmin-4, ymax+10);
    MoveTo(xmin, ymax); LineTo(xmin+4, ymax+10);
    MoveTo(xmin, ymin); LineTo(xmin, ymax+9);

    MoveTo(xmax, ymin); LineTo(xmax-10, ymin+4);
    MoveTo(xmax, ymin); LineTo(xmax-10, ymin-4);
    MoveTo(xmin, ymin); LineTo(xmax-9, ymin);
    Canvas.TextOut(xmin+8, ymax, ytitle);
    Canvas.TextOut(xmax-5*Length(xtitle)-3, ymin+Font.Height-8, xtitle);

    guess_grid(xvalmin, xvalmax, step, xvmin, xvmax);
    v := xvmin + step;
    while v < xvmax do
    begin
      c := round(xmin + (v-xvmin)*(xmax-xmin)/(xvmax-xvmin));
      s := Format('%.4g', [v]);
      TextOut(c-3*Length(s), ymin+Font.Height-8, s);
      MoveTo(c, ymin-3); LineTo(c, ymin+3);
      v := v + step;
    end;

    guess_grid(yvalmin, yvalmax, step, yvmin, yvmax);
    v := yvmin + step;
    while v < yvmax do
    begin
      c := round(ymin + (v-yvmin)*(ymax-ymin)/(yvmax-yvmin));
      s := Format('%.4g', [v]);
      TextOut(xmin+8, c-5, s);
      MoveTo(xmin-3, c); LineTo(xmin+3, c);
      v := v + step;
    end;

    xmin := round(xmin + (xvalmin-xvmin)*(xmax-xmin)/(xvmax-xvmin));
    xmax := round(xmin + (xvalmax-xvmin)*(xmax-xmin)/(xvmax-xvmin));

    ymin := round(ymin + (yvalmin-yvmin)*(ymax-ymin)/(yvmax-yvmin));
    ymax := round(ymin + (yvalmax-yvmin)*(ymax-ymin)/(yvmax-yvmin));
  end;
end;

procedure DrawGraph(Canvas: TCanvas; YY: Graph_Array;
                    yvmin, yvmax, xvmin, xvmax: double;
                    ytitle, xtitle: string);
var
  i, x, y, ymin, ymax, xmin, xmax: integer;
  xscale, yscale: double;
begin
  with Canvas do
  begin
    ymin := ClipRect.Bottom - 10;
    ymax := ClipRect.Top + 10;
    xmin := ClipRect.Left + 10;
    xmax := ClipRect.Right - 10;
    Brush.Color := RGB(255,255,255);
    Pen.Color := RGB(0,0,0);
    FillRect(ClipRect);
    DrawGrid(Canvas, xtitle, ytitle, xmin, ymin, xmax, ymax, xvmin, yvmin, xvmax, yvmax);

    xscale := (xmax-xmin)/(Graph_N-0.5);
    yscale := (ymax-ymin)/(yvmax-yvmin);

    Brush.Color := RGB(0,0,0);
    Pen.Color := RGB(180,0,0);
    y := round(YY[0] * yscale + ymin);
    x := round(0.5 * xscale + xmin);
    MoveTo(round(xmin), round(ymin));
    LineTo(round(xmin), y);
    LineTo(x, y);
    for i := 1 to GR_N-1 do
    begin
      x := round((i-0.5) * xscale + xmin);
      y := round(YY[i] * yscale + ymin);
      LineTo(x, y);
      x := round((i+0.5) * xscale + xmin);
      LineTo(x, y);
    end;
    LineTo(round(xmax), round(ymin));
  end;
end;


{ ----------------------- User Interface code ----------------------- }

procedure TRunner.Execute;
var time_step: double;
    T, dT: TDateTime;
const frame_dt = 1/(24*60*60*40);
begin
  while not Form.StopFlag do
  begin
    time_step := Form.time_step;
    T := Time;
    dT := 0;
    while dT < frame_dt do
    begin
      make_step(25, time_step);
      dT := Time-T;
    end;
    Synchronize(Form.RedrawImage);
  end;
end;

procedure TMainForm.FormCreate(Sender: TObject);
begin
  //randomize;
  mode := Mode_Fixed;
  initial_K := initial_kinetic;
  reset_model(initial_x1, initial_x2, initial_y1, initial_y2, initial_K);

  time_step := initial_dt;
  DecimalSeparator := '.';
  Running := False;
  ParticlesView.Caption := IntToStr(N);
  TimeStepEdit.Text := FloatToStr(time_step);
  SPF := 100;
  SPFEdit.Text := IntToStr(SPF);

  Potential_EEdit.Caption := FloatToStr(potential_e);
  Potential_DEdit.Caption := FloatToStr(potential_d);
  CellsizeXEdit.Caption := FloatToStr(cellsize_x);
  CellsizeYEdit.Caption := FloatToStr(cellsize_y);
  InitialKEdit.Text := FloatToStr(initial_K);
  self.RedrawImage;
end;

procedure TMainForm.RedrawImage;
var ymin, ymax, xmin, xmax: double;
begin
  DrawAtoms(AtomView.Canvas);
  if AdditionalPages.ActivePageIndex = 1 then
  begin
    UpdateStatistics;
  end else begin
    case AdditionalPages.ActivePageIndex of
      2: begin
           calc_GR(GR, ymin, ymax, xmin, xmax);
           DrawGraph(GRView.Canvas, GR, ymin, ymax, xmin, xmax,
                     'g(r)', 'r/d');
          end;
      3: begin
           calc_nK(GR, ymin, ymax, xmin, xmax);
           DrawGraph(KEnergyView.Canvas, GR, ymin, ymax, xmin, xmax,
                     'P, %', 'K, K');
          end;
      4: begin
           calc_nU(GR, ymin, ymax, xmin, xmax);
           DrawGraph(UEnergyView.Canvas, GR, ymin, ymax, xmin, xmax,
                     'P, %', 'U, K');
          end;
      5: begin
           calc_nE(GR, ymin, ymax, xmin, xmax);
           DrawGraph(EnergyView.Canvas, GR, ymin, ymax, xmin, xmax,
                     'P, %', 'E, K');
          end;
    end;
  end;
end;

procedure TMainForm.DrawAtoms(Canvas: TCanvas);
var
  i, x, y: integer;
  xscale, xshift, yscale, yshift: double;
begin
  with Canvas do
  begin
    xscale := (ClipRect.Right-ClipRect.Left-20)/cellsize_x;
    yscale := (ClipRect.Bottom-ClipRect.Top-20)/cellsize_y;
    xshift := ClipRect.Left+10;
    yshift := ClipRect.Top+10;
    Brush.Color := RGB(255,255,255);
    Pen.Color := RGB(140,140,140);
    FillRect(ClipRect);
    MoveTo(ClipRect.Left+10, ClipRect.Top+10);
    LineTo(ClipRect.Right-10, ClipRect.Top+10);
    LineTo(ClipRect.Right-10, ClipRect.Bottom-10);
    LineTo(ClipRect.Left+10, ClipRect.Bottom-10);
    LineTo(ClipRect.Left+10, ClipRect.Top+10);
    Brush.Color := RGB(255,255,255);
    Pen.Color := RGB(155,155,155);
    for i:= 1 to N do
    begin
      x := round(atoms[i].x * xscale + xshift);
      y := round(atoms[i].y * yscale + yshift);
      RoundRect(x-2, y-2, x+3, y+3, 2, 2);
    end;
  end;
end;

procedure TMainForm.UpdateStatistics;
var kmin, kmax, kave, umin, umax, uave, emin, emax, eave, etotal: double;
    dr2, vx2, vy2: double;
begin
  calc_energy(atoms, EA, kmin, kmax, kave, umin, umax, uave, emin, emax, eave, etotal);
  calc_deltas(dr2, vx2, vy2);

  TotalTimeView.Caption := FloatToStrF(total_time, ffGeneral, 10, 4);
  TotalEnergyView.Caption := FloatToStrF(etotal, ffGeneral, 10, 4);
  MinKEnergyView.Caption := FloatToStrF(kmin, ffGeneral, 6, 4);
  AveKEnergyView.Caption := FloatToStrF(kave, ffGeneral, 6, 4);
  MaxKEnergyView.Caption := FloatToStrF(kmax, ffGeneral, 6, 4);
  MinUEnergyView.Caption := FloatToStrF(umin, ffGeneral, 6, 4);
  AveUEnergyView.Caption := FloatToStrF(uave, ffGeneral, 6, 4);
  MaxUEnergyView.Caption := FloatToStrF(umax, ffGeneral, 6, 4);
  ErrorView.Caption := FloatToStrF(abs((etotal-initial_energy)/initial_energy), ffExponent, 1, 1);
  DR2View.Caption := FloatToStrF(dr2, ffGeneral, 6, 4);
  if total_time > 0 then
  begin
    DView.Caption := FloatToStrF(dr2/2.0/total_time, ffGeneral, 6, 4);
    XPressureView.Caption := FloatToStrF(impulse_x/cellsize_y/total_time, ffGeneral, 6, 4);
    YPressureView.Caption := FloatToStrF(impulse_y/cellsize_x/total_time, ffGeneral, 6, 4);
  end else begin
    DView.Caption := '0/0';
    XPressureView.Caption := '0/0';
    YPressureView.Caption := '0/0';
  end;
end;

procedure TMainForm.StopRunning;
begin
  if Running then
  begin
    StopFlag := True;
    Runner.Priority := tpNormal;
    Runner.WaitFor;
    Runner.Destroy;
    Runner := nil;
    Running := False;
    StartButton.Caption := 'Go';
  end;
end;

procedure TMainForm.StartRunning;
begin
  if not Running then
  begin
    StopFlag := False;
    Runner := TRunner.Create(true);
    (Runner as TRunner).Form := self;
    Runner.Priority := tpLower;
    Running := True;
    Runner.Resume;
    StartButton.Caption := 'Stop';
  end;
end;

procedure TMainForm.StartButtonClick(Sender: TObject);
begin
  if not Running then
  begin
    StartButton.Enabled := False;
    StartRunning;
    StartButton.Enabled := True;
  end else begin
    StartButton.Enabled := False;
    StopRunning;
    StartButton.Enabled := True;
  end;
end;

procedure TMainForm.FormDestroy(Sender: TObject);
begin
  StopRunning;
end;

procedure TMainForm.ResetButtonClick(Sender: TObject);
begin
  if Running then
  begin
    StopRunning;
    reset_model(initial_x1, initial_x2, initial_y1, initial_y2, initial_K);
    self.RedrawImage;
    StartRunning;
  end else begin
    reset_model(initial_x1, initial_x2, initial_y1, initial_y2, initial_K);
    self.RedrawImage;
  end;
end;

procedure TMainForm.QuitButtonClick(Sender: TObject);
begin
  Application.Terminate;
end;

procedure TMainForm.FixedBorderSwitchClick(Sender: TObject);
begin
  if Running then
  begin
    StopRunning;
    mode := Mode_Fixed;
    StartRunning;
  end else mode := Mode_Fixed;
end;

procedure TMainForm.PeriodicBorderSwitchClick(Sender: TObject);
begin
  if Running then
  begin
    StopRunning;
    mode := Mode_Periodic;
    StartRunning;
  end else mode := Mode_Periodic;
end;

procedure TMainForm.SPFApplyButtonClick(Sender: TObject);
var steps: integer;
begin
  try
    steps := StrToInt(SPFEdit.Text);
    SPF := steps;
  except
    on EConvertError do
    begin
      Application.MessageBox('That is not a number', 'Validation error');
    end;
  end;
end;

procedure TMainForm.OneFrameButtonClick(Sender: TObject);
begin
  if not Running then
  begin
    OneFrameButton.Enabled := False;
    make_step(SPF, time_step);
    OneFrameButton.Enabled := True;
    RedrawImage;
  end else begin
    StopRunning;
    OneFrameButtonClick(Sender);
  end;
end;

procedure TMainForm.AdditionalPagesChange(Sender: TObject);
begin
  RedrawImage;
end;

procedure TMainForm.ApplyTimeStepButtonClick(Sender: TObject);
var dt: double;
begin
  try
    dt := StrToFloat(TimeStepEdit.Text);
    time_step := dt;
  except
    on EConvertError do
    begin
      Application.MessageBox('That is not a number', 'Validation error');
    end;
  end;
end;

procedure TMainForm.InitialKApplyButtonClick(Sender: TObject);
var K: double;
begin
  try
    K := StrToFloat(InitialKEdit.Text);
    initial_K := K;
  except
    on EConvertError do
    begin
      Application.MessageBox('That is not a number', 'Validation error');
    end;
  end;
end;


end.
