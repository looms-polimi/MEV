package MEV "Models of the Multiple Emergency Ventilator"
  extends Modelica.Icons.Package;

  package Info "Information"
    extends Modelica.Icons.Information;
    annotation(
      DocumentationClass = true,
      Documentation(info = "<html><head></head><body><p>
This package was created from scratch in the last week of March 2020 to support
the system and control design of the Multiple Emergency Ventilator, a 10+ bed
low-cost medical ventilatrion system to address the needs of pandemic outbreaks
like COVID-19, that require a sudden expansion of intensive care units and
of mechanical ventilation of sedated and curarized patients.</p>
<p>The aim of the package is twofold:
</p><ul>
<li>To assist in the design of the control strategy for the MEV system.</li>
<li>To demonstrate how the Modelica language, an open standard developed by the
non-profit Modelica Association, can be effective in the fast prototyping of
medical devices.</li>
</ul><p></p>
<p>The models contained in this package were run successfully tested with OpenModelica 1.16.0 and with Dymola 2020x, obtaining the same results.</p>
<p>For further information about the MEV project, please refer to the
<a href=\"https://mev.deib.polimi.it/\">mev.deib.polimi.it</a> website.</p>
<p>Main author: Francesco Casella, <a href=\"mailto:francesco.casella@polimi.it\">
francesco.casella@polimi.it</a>.
</p></body></html>"));
  end Info;

  model Ambient
    constant Modelica.SIunits.SpecificHeatCapacity Rstar = Modelica.Constants.R / MM "Gas constant for air";
    constant Modelica.SIunits.MolarMass MM = 0.029 "Molar mass of air";
    constant Modelica.SIunits.PerUnit gamma = 1.4 "cp/cv";
    parameter Types.AbsolutePressure p = 101325 "Ambient absolute pressure";
    parameter Modelica.SIunits.Temperature T = 298.15 "Reference gas temperature";
    parameter Boolean useOnOffControl = false "=true uses on-off controller, otherwise linear controller";
    annotation(
      Icon(coordinateSystem(preserveAspectRatio = false), graphics = {Rectangle(extent = {{-100, 100}, {100, -100}}, lineColor = {0, 0, 0}), Text(extent = {{-80, 64}, {80, -58}}, lineColor = {28, 108, 200}, textString = "Ambient")}),
      Diagram(coordinateSystem(preserveAspectRatio = false)),
      defaultComponentName = "ambient",
      defaultComponentPrefixes = "inner",
      missingInnerMessage = "The Ambient object is missing");
  end Ambient;

  package Interfaces "Physical connectors"
    extends Modelica.Icons.InterfacesPackage;

    connector PneumaticPort
      flow Types.MassFlowRate w "Mass flow rate flowing into port";
      Types.RelativePressure p "Relative Pressure w.r.t. ambient";
    end PneumaticPort;

    connector PneumaticPort_a
      extends PneumaticPort;
      annotation(
        Icon(graphics = {Ellipse(lineColor = {111, 164, 171}, fillColor = {166, 244, 255}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}, endAngle = 360)}));
    end PneumaticPort_a;

    connector PneumaticPort_b
      extends PneumaticPort;
      annotation(
        Icon(graphics = {Ellipse(lineColor = {111, 164, 171}, fillColor = {166, 244, 255}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}, endAngle = 360), Ellipse(extent = {{-80, 80}, {80, -80}}, lineColor = {111, 164, 171}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid)}));
    end PneumaticPort_b;
  end Interfaces;

  package Sources "Source components"
    extends Modelica.Icons.Package;

    model PressureSource
      outer Ambient ambient "Ambient model";
      parameter Types.AbsolutePressure pabs = ambient.p "Absolute pressure" annotation(
        Dialog(enable = not useRelativePressure));
      parameter Types.RelativePressure prel(displayUnit = "Pa") = 0 "Relative pressure" annotation(
        Dialog(enable = useRelativePressure));
      parameter Boolean useRelativePressure = true "= true if relative pressure is to be used";
      Interfaces.PneumaticPort_a port annotation(
        Placement(visible = true, transformation(origin = {-100, 4}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      if useRelativePressure then
        port.p = prel;
      else
        port.p = pabs - ambient.p;
      end if;
      annotation(
        Icon(graphics = {Ellipse(origin = {-1, 0}, lineColor = {131, 189, 197}, fillColor = {166, 244, 255}, fillPattern = FillPattern.Solid, extent = {{-99, 100}, {101, -100}}, endAngle = 360), Text(extent = {{-140, -104}, {140, -140}}, lineColor = {28, 108, 200}, textString = "%name")}));
    end PressureSource;

    model FlowSource
      outer Ambient ambient "Ambient model";
      Interfaces.PneumaticPort_b port annotation(
        Placement(visible = true, transformation(origin = {100, -4}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput q "Volume flow rate in l/min" annotation(
        Placement(visible = true, transformation(origin = {-112, -6}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-100, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      Media.Air air "Air model";
      Types.MassFlowRate w "Leaving mass flow rate";
    equation
      air.p = port.p + ambient.p;
      w = air.rho * q / (1000 * 60);
      port.w = -w;
      annotation(
        Icon(graphics = {Rectangle(origin = {1, 0}, lineColor = {131, 189, 197}, fillColor = {166, 244, 255}, fillPattern = FillPattern.Solid, extent = {{-101, 40}, {99, -40}}), Polygon(origin = {-4, -8}, fillColor = {255, 255, 255}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, points = {{-16, 30}, {-16, -18}, {44, 6}, {-16, 30}}), Text(extent = {{-138, -46}, {140, -82}}, lineColor = {28, 108, 200}, textString = "%name")}));
    end FlowSource;
  end Sources;

  package Components "Basic components"
    extends Modelica.Icons.Package;

    model Cylinder
      Interfaces.PneumaticPort_a inlet annotation(
        Placement(visible = true, transformation(origin = {-102, -98}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-60, -80}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      Interfaces.PneumaticPort_b outlet annotation(
        Placement(visible = true, transformation(origin = {100, -98}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {60, -80}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      outer Ambient ambient;
      parameter Modelica.SIunits.Area A "Cylinder cross section";
      parameter Modelica.SIunits.Volume Vref "Cylinder volume at reference level";
      parameter Modelica.SIunits.Length yref "Reference level";
      parameter Modelica.SIunits.Length ystart "Start value of the level";
      parameter Modelica.SIunits.Length ymin = 1e-3 "Minimum level";
      parameter Modelica.SIunits.Length ymax "Maximum level";
      parameter Types.RelativePressure pstart = 0 "Start value of relative pressure";
      Modelica.SIunits.Mass M(stateSelect = StateSelect.prefer) "Mass";
      Modelica.SIunits.Length y(stateSelect = StateSelect.prefer) "Level";
      Types.Volume V "Volume";
      Types.RelativePressure p "Relative air pressure in the cylinder";
      Modelica.SIunits.Force F "Net force applied to the piston, discounting external atmospheric pressure";
      Types.MassFlowRate w_in "Entering mass flow rate";
      Types.MassFlowRate w_out "Leaving mass flow rate";
      final Types.AbsolutePressure pamb = ambient.p "Ambient pressure";
      Media.AirAdiabatic air;
      parameter Boolean init_y = false "Initialize the level";
      parameter Boolean init_p = false "Initialize the pressure";
      Modelica.Mechanics.Translational.Interfaces.Flange_a pistonFlange annotation(
        Placement(transformation(extent = {{-10, -10}, {10, 10}}), iconTransformation(extent = {{-10, 10}, {10, 30}})));
    initial equation
      if init_y then
        y = ystart;
      end if;
      if init_p then
        p = pstart;
      end if;
    equation
      der(M) = w_in - w_out;
      V = A * y;
      M = air.rho * V;
      air.p = p + pamb;
      F = p * A;
// Boundary condition
      inlet.p = p;
      outlet.p = p;
      w_in = inlet.w;
      w_out = -outlet.w;
      pistonFlange.f = -F;
      pistonFlange.s = y;
// Assertions
      assert(y > ymin, "The level is below the minimum value");
      assert(y < ymax, "The level is above the maximum value");
      annotation(
        Icon(graphics = {Rectangle(origin = {0.6, -45}, fillColor = {166, 244, 255}, fillPattern = FillPattern.Solid, extent = {{-60.6, 55}, {59.4, -55}}, lineColor = {0, 0, 0}), Rectangle(origin = {0, -10}, lineThickness = 1, extent = {{-60, 90}, {60, -90}}), Rectangle(extent = {{-58, 10}, {58, 10}}, lineColor = {0, 0, 0}, fillColor = {136, 136, 136}, fillPattern = FillPattern.Solid, lineThickness = 0.5)}));
    end Cylinder;

    model Pipe
      Interfaces.PneumaticPort_a inlet annotation(
        Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Interfaces.PneumaticPort_b outlet annotation(
        Placement(visible = true, transformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      parameter Modelica.SIunits.Length D "Pipe diameter";
      parameter Modelica.SIunits.Length L "Pipe length";
      parameter Modelica.SIunits.PerUnit cf = 0.005 "Fanning friction factor";
      parameter Types.Velocity unom "Nominal air speed";
      parameter Modelica.SIunits.PerUnit delta = 0.1 "Fraction of nominal velocity for friction regularization";
      parameter Types.RelativePressure pstart = 0 "start value for relative pressure p";
      constant Modelica.SIunits.PerUnit pi = Modelica.Constants.pi;
      final parameter Modelica.SIunits.Area A = pi * D ^ 2 / 4 "Pipe cross-sectional surface";
      final parameter Modelica.SIunits.Volume V = A * L;
      final parameter Modelica.SIunits.Length omega = pi * D "Wet perimeter";
      final parameter Modelica.SIunits.Length l = L / 2 "Length of a single finite vol.";
      Media.AirAdiabatic air;
      Types.MassFlowRate win "Entering flow rate";
      Types.MassFlowRate wout "Leaving flow rate";
      Types.Velocity uin "Entering fluid speed";
      Types.Velocity uout "Leaving fluid velocity";
      Types.RelativePressure pin "Entering relative pressure";
      Types.RelativePressure pout "Leaving relative pressure";
      Types.RelativePressure dpin "Pressure loss, inlet side";
      Types.RelativePressure dpout "Pressure loss, outlet side";
      Types.RelativePressure dp "Pressure loss, total";
      Types.VelocitySquared uin2 "Regularized squared inlet velocity";
      Types.VelocitySquared uout2 "Regularized squared inlet velocity";
      Types.RelativePressure p(stateSelect = StateSelect.prefer) "Relative pressure at the middle of the pipe";
      Modelica.SIunits.Mass M(stateSelect = StateSelect.avoid) "Mass of fluid in the pipe";
    protected
      function regSquare = Modelica.Fluid.Utilities.regSquare;
    equation
// Mass balance
      der(M) = win - wout;
      M = air.rho * V;
// Momentum balance
      uin = win / (air.rho * A);
      uout = wout / (air.rho * A);
      uin2 = homotopy(regSquare(uin, unom * delta), uin * unom);
      uout2 = homotopy(regSquare(uout, unom * delta), uout * unom);
      dpin = pin - p;
      dpout = p - pout;
      dp = dpin + dpout;
      A*dpin = cf/2*air.rho*omega*uin2*l "Momentum balance, entering side";
      A*dpout= cf/2*air.rho*omega*uout2*l "Momentum balance, leaving side";
// Air model
      air.p = p;
// Boundary values
      win = inlet.w;
      wout = -outlet.w;
      pin = inlet.p;
      pout = outlet.p;
    initial equation
      p = pstart;
      annotation(
        Icon(graphics = {Rectangle(origin = {0, 3}, fillColor = {166, 244, 255}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{-100, 17}, {100, -23}}, lineColor = {131, 189, 197})}));
    end Pipe;

    model Compliance
      outer Ambient ambient "Ambient model";
      parameter Types.Compliance C "Capacitance";
      parameter Types.RelativePressure pstart = 0 "Initial relative pressure";
      Interfaces.PneumaticPort_a port annotation(
        Placement(visible = true, transformation(origin = {-100, 4}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Types.RelativePressure p "Relative pressure";
      Types.MassFlowRate w "Absolute pressure";
      Types.Volume V "Volume of breathed air";
      Types.Volume Vb "Total volume of breathed air";
      Media.Air air;
    initial equation
      p = pstart;
      V = 0;
      Vb = 0;
    equation
      C * der(p) = w / air.rho;
      der(V) = w / air.rho;
      der(Vb) = max(der(V), 0);
// Reference air properties
      air.p = ambient.p;
// Boundary conditions
      w = port.w;
      p = port.p;
      annotation(
        Icon(graphics = {Ellipse(origin = {-1, 0}, lineColor = {131, 189, 197}, fillColor = {166, 244, 255}, fillPattern = FillPattern.Solid, extent = {{-99, 100}, {101, -100}}, endAngle = 360), Ellipse(extent = {{-60, 80}, {60, -80}}, lineColor = {131, 189, 197}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid)}));
    end Compliance;

    model Resistance
      Interfaces.PneumaticPort_a inlet annotation(
        Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Interfaces.PneumaticPort_b outlet annotation(
        Placement(visible = true, transformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      outer Ambient ambient "Ambient model";
      parameter Types.Resistance R "Resistace";
      Media.AirAdiabatic air;
      Types.MassFlowRate w "Flow rate";
      Types.RelativePressure pin "Entering relative pressure";
      Types.RelativePressure pout "Leaving relative pressure";
      Types.RelativePressure dp "Pressure loss, inlet side";
      Types.VolumeFlowRate q "Volume flow rate";
    protected
      function regSquare = Modelica.Fluid.Utilities.regSquare;
    equation
      dp = pin - pout;
      q = w / air.rho;
      dp = R * q;
// Air model
      air.p = pin + ambient.p;
// Boundary values
      w = inlet.w;
      w = -outlet.w;
      pin = inlet.p;
      pout = outlet.p;
      annotation(
        Icon(graphics = {Line(points = {{-60, 0}, {-100, 0}}, color = {0, 0, 0}, thickness = 1), Line(points = {{100, 0}, {60, 0}}, color = {0, 0, 0}, thickness = 1), Rectangle(extent = {{-60, 20}, {62, -20}}, lineColor = {0, 0, 0}, lineThickness = 1)}));
    end Resistance;

    model PressureLoss
      Interfaces.PneumaticPort_a inlet annotation(
        Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Interfaces.PneumaticPort_b outlet annotation(
        Placement(visible = true, transformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      outer Ambient ambient "Ambient model";
      parameter Types.MassFlowRate wnom "Nominal mass flow rate";
      parameter Types.RelativePressure dpnom "Nominal pressure drop";
      parameter Types.Density rhonom(displayUnit = "kg/m3") = 1.2 "Nominal air density";
      parameter Modelica.SIunits.PerUnit delta = 0.1 "Fraction of nominal mass flow rate for friction regularization";
      Media.AirAdiabatic air;
      Types.MassFlowRate w "Flow rate";
      Types.RelativePressure pin "Entering relative pressure";
      Types.RelativePressure pout "Leaving relative pressure";
      Types.RelativePressure dp "Pressure loss, inlet side";
    protected
      function regSquare = Modelica.Fluid.Utilities.regSquare;
    equation
      dp = pin - pout;
      dp = homotopy(regSquare(w, wnom * delta), dpnom / wnom * w);
// Air model
      air.p = pin + ambient.p;
// Boundary values
      w = inlet.w;
      w = -outlet.w;
      pin = inlet.p;
      pout = outlet.p;
      annotation(
        Icon(graphics = {Polygon(lineThickness = 1, points = {{-100, 100}, {-100, -100}, {100, 100}, {100, -100}, {-100, 100}})}));
    end PressureLoss;

    model ControlValve
      Interfaces.PneumaticPort_a inlet annotation(
        Placement(visible = true, transformation(origin = {-102, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      Interfaces.PneumaticPort_b outlet annotation(
        Placement(visible = true, transformation(origin = {100, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput opening(unit = "1") annotation(
        Placement(visible = true, transformation(origin = {-4, 64}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {0, 42}, extent = {{-20, -20}, {20, 20}}, rotation = -90)));
      outer Ambient ambient "Ambient model";
      parameter Modelica.SIunits.Area Av "Maximum flow coefficient";
      parameter Types.MassFlowRate wnom "Nominal mass flow rate";
      parameter Types.RelativePressure dpnom "Nominal delta-p";
      parameter Modelica.SIunits.PerUnit delta = 0.01 "Fraction of nominal delta-p for sqrt regularization";
      final parameter Real Cv = Av * 1e6 / 24 "Cv coefficient of the valve";
      final parameter Real Kv = Av * 1e6 / 27.7 "Kv coefficient of the valve";
      Media.Air air "inlet air properties";
      Types.MassFlowRate w "Mass flow rate";
      Types.RelativePressure pin "Entering pressure, relative";
      Modelica.SIunits.Pressure pout "Leaving pressure, relative";
      Types.RelativePressure dp "Valve delta-p";
      Types.VolumeFlowRate q "Volume flow rate";
    protected
      function regRoot = Modelica.Fluid.Utilities.regRoot;
    equation
      dp = pin - pout;
      w = homotopy(Av * opening * sqrt(air.rho) * regRoot(dp, dpnom * delta), wnom / dpnom * dp * opening);
      q = w / air.rho;
      air.p = pin + ambient.p;
// Boundary conditions
      inlet.w = w;
      outlet.w = -w;
      inlet.p = pin;
      outlet.p = pout;
      annotation(
        Icon(graphics = {Polygon(lineThickness = 1, points = {{-100, 100}, {-100, -100}, {100, 100}, {100, -100}, {-100, 100}}), Text(extent = {{-140, -108}, {140, -146}}, lineColor = {28, 108, 200}, textString = "%name")}));
    end ControlValve;
  end Components;

  package ParameterizedComponents "Sized components for the MEV system"
    extends Modelica.Icons.Package;

    model BellJar
      outer Ambient ambient "Ambient conditions";
      constant Modelica.SIunits.Acceleration g = Modelica.Constants.g_n "Acceleration of gravity";
      parameter Modelica.SIunits.Area A = 20e-3 / 0.174 "Cylinder cross section";
      parameter Modelica.SIunits.Length ystart = (ymin + ymax) / 2 "Start value of the level";
      parameter Modelica.SIunits.Length ymax = 0.174 "Maximum level";
      parameter Modelica.SIunits.Length ymin = 0 "Minimum level";
      parameter Modelica.SIunits.Length yref = ymin "Reference level";
      parameter Modelica.SIunits.Volume Vref = 15.27e-3 "Volume at reference level";
      parameter Types.RelativePressure pdes = 2400 "Relative supply pressure design value";
      final parameter Modelica.SIunits.Mass M = pdes * A / g "Piston weight";
      final parameter Types.AbsolutePressure pa = ambient.p "Absolute ambient pressure";
      parameter Modelica.SIunits.TranslationalDampingConstant d = 10 "Damping constant";
      Components.Cylinder cylinder(A = A, Vref = Vref, yref = yref, ystart = ystart, ymax = ymax, ymin = ymin, pstart = pstart, init_y = true, init_p = true) annotation(
        Placement(transformation(extent = {{-40, -62}, {0, -22}})));
      Modelica.Mechanics.Translational.Components.Mass piston(m = M, v(fixed = true)) annotation(
        Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = -90, origin = {-20, -4})));
      Modelica.Mechanics.Translational.Sources.Force force annotation(
        Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = -90, origin = {-20, 26})));
      Modelica.Mechanics.Translational.Components.Fixed fixed annotation(
        Placement(transformation(extent = {{10, -44}, {30, -24}})));
      Modelica.Mechanics.Translational.Components.Damper damper(d = d) annotation(
        Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = 90, origin = {20, -12})));
      Modelica.Blocks.Sources.RealExpression forceValue(y = -M * g) annotation(
        Placement(transformation(extent = {{-66, 40}, {-36, 60}})));
      Interfaces.PneumaticPort_a inlet annotation(
        Placement(visible = true, transformation(origin = {-56, -58}, extent = {{-6, -6}, {6, 6}}, rotation = 0), iconTransformation(origin = {-60, -80}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      Interfaces.PneumaticPort_b outlet annotation(
        Placement(visible = true, transformation(origin = {14, -58}, extent = {{-6, -6}, {6, 6}}, rotation = 0), iconTransformation(origin = {60, -80}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      Modelica.Mechanics.Translational.Sensors.PositionSensor heightSensor annotation(
        Placement(transformation(extent = {{34, -4}, {54, 16}})));
      Modelica.Blocks.Interfaces.RealOutput y "Height of piston [m]" annotation(
        Placement(transformation(extent = {{64, -4}, {84, 16}}), iconTransformation(extent = {{-72, 50}, {-92, 70}})));
      parameter Types.RelativePressure pstart = pdes "Start value of relative pressure";
    equation
      connect(piston.flange_b, cylinder.pistonFlange) annotation(
        Line(points = {{-20, -14}, {-20, -38}}, color = {0, 127, 0}));
      connect(force.flange, piston.flange_a) annotation(
        Line(points = {{-20, 16}, {-20, 6}}, color = {0, 127, 0}));
      connect(fixed.flange, damper.flange_a) annotation(
        Line(points = {{20, -34}, {20, -22}}, color = {0, 127, 0}));
      connect(damper.flange_b, piston.flange_a) annotation(
        Line(points = {{20, -2}, {20, 6}, {-20, 6}}, color = {0, 127, 0}));
      connect(forceValue.y, force.f) annotation(
        Line(points = {{-34.5, 50}, {-20, 50}, {-20, 38}}, color = {0, 0, 127}));
      connect(inlet, cylinder.inlet) annotation(
        Line(points = {{-56, -58}, {-32, -58}}, color = {111, 164, 171}));
      connect(inlet, inlet) annotation(
        Line(points = {{-56, -58}, {-56, -58}}, color = {111, 164, 171}));
      connect(cylinder.outlet, outlet) annotation(
        Line(points = {{-8, -58}, {14, -58}}, color = {111, 164, 171}));
      connect(heightSensor.s, y) annotation(
        Line(points = {{55, 6}, {74, 6}}, color = {0, 0, 127}));
      connect(piston.flange_a, heightSensor.flange) annotation(
        Line(points = {{-20, 6}, {34, 6}}, color = {0, 127, 0}));
      annotation(
        Icon(coordinateSystem(extent = {{-80, -100}, {80, 100}}), graphics = {Rectangle(origin = {0.6, -45}, fillColor = {166, 244, 255}, fillPattern = FillPattern.Solid, extent = {{-60.6, 55}, {59.4, -55}}, lineColor = {0, 0, 0}), Rectangle(origin = {0, -10}, lineThickness = 1, extent = {{-60, 90}, {60, -90}}), Rectangle(extent = {{-60, 32}, {60, 8}}, lineColor = {0, 0, 0}, fillColor = {136, 136, 136}, fillPattern = FillPattern.Solid, lineThickness = 0.5), Line(points = {{-40, 86}, {-40, 32}, {-40, 34}}, color = {0, 0, 0}), Line(points = {{-44, 40}, {-36, 40}}, color = {0, 0, 0}), Line(points = {{-44, 46}, {-36, 46}}, color = {0, 0, 0}), Line(points = {{-44, 52}, {-36, 52}}, color = {0, 0, 0}), Line(points = {{-44, 58}, {-36, 58}}, color = {0, 0, 0}), Line(points = {{-44, 64}, {-36, 64}}, color = {0, 0, 0}), Line(points = {{-44, 70}, {-36, 70}}, color = {0, 0, 0}), Line(points = {{-72, 60}, {-40, 60}, {-40, 60}}, color = {0, 0, 0})}),
        Diagram(coordinateSystem(extent = {{-80, -100}, {80, 100}})));
    end BellJar;

    model PipeSegment
      extends Components.Pipe(D = 0.050, cf = 0.01, unom = 1, pstart = 2400, L = 2.2);
    end PipeSegment;

    model PatientValve
      extends Components.ControlValve(dpnom = 100, wnom = 0.3e-3, Av = 0.5 * 1e-4);
    end PatientValve;

    model SupplyValve
      extends Components.ControlValve(dpnom = 100, wnom = 0.3e-3, Av = 5.75e-6);
    end SupplyValve;

    model StandardPatient
      PatientValve valve_in annotation(
        Placement(visible = true, transformation(origin = {-60, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Components.Resistance resistance(R = kR * R) annotation(
        Placement(transformation(extent = {{-20, -20}, {20, 20}}, rotation = -90, origin = {12, -54})));
      Sources.PressureSource peep(prel = 1200) annotation(
        Placement(visible = true, transformation(origin = {88, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      PatientValve valve_out annotation(
        Placement(visible = true, transformation(origin = {42, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Components.Compliance compliance(C = kC * C, pstart = pstart) annotation(
        Placement(transformation(extent = {{-8, -128}, {32, -88}})));
      Interfaces.PneumaticPort_a supply annotation(
        Placement(transformation(extent = {{-120, 10}, {-100, 30}}), iconTransformation(extent = {{-20, 60}, {20, 100}})));
      Modelica.Blocks.Math.Add add(k2 = -1) annotation(
        Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = -90, origin = {42, 50})));
      Modelica.Blocks.Sources.Constant one(k = 1) annotation(
        Placement(transformation(extent = {{84, 70}, {64, 90}})));
      Modelica.Blocks.Interfaces.RealInput valveCommand annotation(
        Placement(transformation(extent = {{-134, 60}, {-94, 100}}), iconTransformation(extent = {{-100, -16}, {-68, 16}})));
      parameter Types.Compliance C = 40/1e6 / (0.01*9.81*999) "Compliance";
      parameter Types.Resistance R = 12/100 * 9.81*999/0.001 "Resistance";
      parameter Types.RelativePressure pstart = 1200 "Initial relative pressure";
      parameter Modelica.SIunits.PerUnit kC = 1 "Correction factor for compliance";
      parameter Modelica.SIunits.PerUnit kR = 1 "Correction factor for resistance";
      Components.ControlValve leak(Av = 2e-4, wnom = 0.01, dpnom = 100) annotation(
        Placement(transformation(extent = {{-10, 10}, {10, -10}}, rotation = -90, origin = {-30, -24})));
      Modelica.Blocks.Sources.RealExpression leakOpening annotation(
        Placement(transformation(extent = {{-88, -38}, {-52, -10}})));
      Sources.PressureSource atmosphere annotation(
        Placement(transformation(extent = {{-40, -72}, {-20, -52}})));
    equation
      connect(add.y, valve_out.opening) annotation(
        Line(points = {{42, 39}, {42, 24.2}}, color = {0, 0, 127}));
      connect(compliance.port, resistance.outlet) annotation(
        Line(points = {{12, -108}, {12, -74}}, color = {111, 164, 171}));
      connect(valve_in.outlet, resistance.inlet) annotation(
        Line(points = {{-50, 20}, {12, 20}, {12, -34}}, color = {111, 164, 171}));
      connect(resistance.inlet, valve_out.inlet) annotation(
        Line(points = {{12, -34}, {12, 20}, {32, 20}}, color = {111, 164, 171}));
      connect(valve_out.outlet, peep.port) annotation(
        Line(points = {{52, 20}, {88, 20}}, color = {111, 164, 171}));
      connect(supply, valve_in.inlet) annotation(
        Line(points = {{-110, 20}, {-70, 20}}, color = {111, 164, 171}));
      connect(add.u1, one.y) annotation(
        Line(points = {{48, 62}, {48, 80}, {63, 80}}, color = {0, 0, 127}));
      connect(valveCommand, valve_in.opening) annotation(
        Line(points = {{-114, 80}, {-60, 80}, {-60, 24.2}}, color = {0, 0, 127}));
      connect(valveCommand, add.u2) annotation(
        Line(points = {{-114, 80}, {36, 80}, {36, 62}}, color = {0, 0, 127}));
      connect(valve_in.outlet, leak.inlet) annotation(
        Line(points = {{-50, 20}, {-30, 20}, {-30, -14}}, color = {111, 164, 171}));
      connect(leak.opening, leakOpening.y) annotation(
        Line(points = {{-34.2, -24}, {-50.2, -24}}, color = {0, 0, 127}));
      connect(leak.outlet, atmosphere.port) annotation(
        Line(points = {{-30, -34}, {-30, -62}}, color = {111, 164, 171}));
      annotation(
        Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-140, -140}, {140, 140}}), graphics = {Line(points = {{0, 60}, {0, 0}, {-10, 0}, {10, 0}, {-10, -14}, {10, -14}, {-10, 0}}, color = {0, 0, 0}), Ellipse(extent = {{-22, -44}, {20, -90}}, lineColor = {0, 0, 0}, fillColor = {255, 255, 255}, fillPattern = FillPattern.None), Line(points = {{0, -14}, {0, -72}, {8, -60}, {40, -60}, {40, -40}}, color = {0, 0, 0}), Line(points = {{40, 34}, {40, -26}, {30, -26}, {50, -26}, {30, -40}, {50, -40}, {30, -26}}, color = {0, 0, 0}), Ellipse(extent = {{26, 60}, {52, 34}}, lineColor = {0, 0, 0}, fillColor = {255, 255, 255}, fillPattern = FillPattern.None)}),
        Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-140, -140}, {140, 140}})));
    end StandardPatient;
  end ParameterizedComponents;

  package SystemModels "Models of the MEV system"
    extends Modelica.Icons.Package;

    model BaseSystem
      parameter Integer N = 10 "Maximum number of patients";
      ParameterizedComponents.BellJar bellJar annotation(
        Placement(transformation(extent = {{-28, -6}, {20, 54}})));
      ParameterizedComponents.SupplyValve supplyValveOnOff annotation(
        Placement(transformation(extent = {{-120, -32}, {-104, -16}})));
      Sources.PressureSource pressureSource(prel = 2e5, useRelativePressure = true) annotation(
        Placement(transformation(extent = {{-178, -10}, {-158, 10}})));
      ParameterizedComponents.PipeSegment pipeSegments[N] annotation(
        Placement(transformation(extent = {{46, -22}, {90, 22}})));
      ParameterizedComponents.StandardPatient patients[N] annotation(
        Placement(transformation(extent = {{78, -80}, {138, -20}})));
      inner Ambient ambient annotation(
        Placement(transformation(extent = {{100, 60}, {120, 80}})));
      ParameterizedComponents.SupplyValve supplyValveModulating annotation(
        Placement(transformation(extent = {{-120, 4}, {-104, 20}})));
      Controls.OnOffControllerWithHysteresis onOffControllerWithHysteresis(level_min = bellJar.ymin + 0.02, level_max = bellJar.ymax - 0.02) annotation(
        Placement(transformation(extent = {{-52, 18}, {-72, 38}})));
      Controls.LinearController linearController(level_sp = (bellJar.ymin + bellJar.ymax) / 2 + 0.030, Kp = 10, T = 1) annotation(
        Placement(transformation(extent = {{-52, 48}, {-72, 68}})));
      MEV.Controls.PatientController patientControllers[N] annotation(
        Placement(visible = true, transformation(origin = {50, -60}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
      Modelica.Blocks.Sources.Step RR[N](each height = 0) "Respiratory rates of patients" annotation(
        Placement(visible = true, transformation(origin = {10, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Step dutyCycle[N](each height = 0) annotation(
        Placement(visible = true, transformation(origin = {10, -76}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(pressureSource.port, supplyValveOnOff.inlet) annotation(
        Line(points = {{-168, 0}, {-140, 0}, {-140, -24}, {-120, -24}}, color = {111, 164, 171}));
      connect(supplyValveOnOff.outlet, bellJar.inlet) annotation(
        Line(points = {{-104, -24}, {-40, -24}, {-40, 0}, {-22, 0}}, color = {111, 164, 171}));
      connect(bellJar.outlet, pipeSegments[1].inlet) annotation(
        Line(points = {{14, 0}, {46, 0}}, color = {111, 164, 171}));
      connect(pipeSegments.outlet, patients.supply) annotation(
        Line(points = {{90, 0}, {108, 0}, {108, -32.8571}}, color = {111, 164, 171}));
      connect(pipeSegments[1:N - 1].outlet, pipeSegments[2:N].inlet) annotation(
        Line(points = {{90, 0}, {46, 0}}, color = {111, 164, 171}));
      connect(pressureSource.port, supplyValveModulating.inlet) annotation(
        Line(points = {{-168, 0}, {-140, 0}, {-140, 12}, {-120, 12}}, color = {111, 164, 171}));
      connect(supplyValveModulating.outlet, bellJar.inlet) annotation(
        Line(points = {{-104, 12}, {-40, 12}, {-40, 0}, {-22, 0}}, color = {111, 164, 171}));
      connect(bellJar.y, onOffControllerWithHysteresis.level) annotation(
        Line(points = {{-28.6, 42}, {-36, 42}, {-36, 28}, {-50, 28}}, color = {0, 0, 127}));
      connect(bellJar.y, linearController.level) annotation(
        Line(points = {{-28.6, 42}, {-36, 42}, {-36, 58}, {-50, 58}}, color = {0, 0, 127}));
      connect(linearController.valveOpening, supplyValveModulating.opening) annotation(
        Line(points = {{-73, 58}, {-112, 58}, {-112, 15.36}}, color = {0, 0, 127}));
      connect(onOffControllerWithHysteresis.valveOpening, supplyValveOnOff.opening) annotation(
        Line(points = {{-73, 28}, {-88, 28}, {-88, -6}, {-112, -6}, {-112, -20.64}}, color = {0, 0, 127}));
  connect(patientControllers.opening, patients.valveCommand) annotation(
        Line(points = {{60, -60}, {75, -60}, {75, -50}, {90, -50}}, color = {0, 0, 127}, thickness = 0.5));
  connect(RR.y, patientControllers.RR) annotation(
        Line(points = {{21, -40}, {28, -40}, {28, -54}, {40, -54}}, color = {0, 0, 127}, thickness = 0.5));
  connect(dutyCycle.y, patientControllers.dutyCycle) annotation(
        Line(points = {{21, -76}, {28, -76}, {28, -66}, {40, -66}}, color = {0, 0, 127}, thickness = 0.5));
      annotation(
        Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-180, -100}, {180, 100}})),
        Diagram(coordinateSystem(extent = {{-200, -100}, {140, 100}}), graphics = {Text(lineColor = {28, 108, 200}, extent = {{52, 20}, {76, 12}}, textString = "%N% X"), Text(lineColor = {28, 108, 200}, extent = {{88, -56}, {112, -64}}, textString = "%N% X")}));
    end BaseSystem;

    model StandardPatientMixSystem "System with standard mix of 10 patients"
      extends MEV.SystemModels.BaseSystem(
        patients(
          C               = {30, 30, 30, 40, 40, 40, 40, 60, 60, 60}/1e6/(0.01*9.81*999),
          R               = { 6, 12, 18,  6, 12, 12, 18,  6, 12, 18}/100*9.81*999/0.001),
        RR(height         = {30, 29, 28, 23, 22, 21, 19, 16, 15, 14},
           each startTime = -10),
        dutyCycle(offset  = {40, 40, 40, 33, 33, 33, 33, 20, 25, 33}));
    end StandardPatientMixSystem;

    model WorstCasePatientSystem "Worst case with 10 patients max C min R"
      extends MEV.SystemModels.BaseSystem(
        patients(
          each C = 60/1e6/(0.01*9.81*999),
          each R = 6/100*9.81*999/0.001),
        RR(each height = 16,
           each startTime = -10),
        dutyCycle(each offset = 20));
    end WorstCasePatientSystem;

    model StandardPatientMix2XSystem "System with 2X standard mix of patients"
      extends MEV.SystemModels.BaseSystem(N = 20,
        patients(
          C = {30, 30, 30, 40, 40, 40, 40, 60, 60, 60,
               30, 30, 30, 40, 40, 40, 40, 60, 60, 60}/1e6/(0.01*9.81*999),
          R = { 6, 12, 18,  6, 12, 12, 18,  6, 12, 18,
                6, 12, 18,  6, 12, 12, 18,  6, 12, 18}/100*9.81*999/0.001),
          RR(height         = {30, 29, 28, 23, 22, 21, 19, 16, 15, 14,
                               30, 29, 28, 23, 22, 21, 19, 16, 15, 14},
             each startTime = -10),
          dutyCycle(offset  = {40, 40, 40, 33, 33, 33, 33, 20, 25, 33,
                               40, 40, 40, 33, 33, 33, 33, 20, 25, 33}));
    end StandardPatientMix2XSystem;
  end SystemModels;

  package Simulations "Simuation scenarios for the MEV system"
    extends Modelica.Icons.ExamplesPackage;

    package LinearControl "Scenarios with linear modulating control of the bell jar"
      extends Modelica.Icons.ExamplesPackage;

      model Scenario1 "System turned on, first patient attached at time = 10, linear control"
        extends Modelica.Icons.Example;
        extends SystemModels.StandardPatientMixSystem(
          RR(startTime = {1e6, 1e6, 1e6, 1e6, 10, 1e6, 1e6, 1e6, 1e6, 1e6}),
          bellJar(ystart = bellJar.ymin + 0.005));
        annotation(
          Diagram(coordinateSystem(extent = {{-200, -100}, {140, 100}})),
          experiment(StopTime = 25, Interval = 0.005, Tolerance = 1e-06, StartTime = 0));
      end Scenario1;

      model Scenario2 "Five patients attached, linear control"
        extends Modelica.Icons.Example;
        extends SystemModels.StandardPatientMixSystem(
          RR(startTime = {-10, 1e6, -10, 1e6, -10, 1e6, 1e6, -10, 1e6, -10}));
        annotation(
          Diagram(coordinateSystem(extent = {{-200, -100}, {140, 100}})),
          experiment(StopTime = 15, Interval = 0.005, Tolerance = 1e-06, StartTime = 0));
      end Scenario2;

      model Scenario3 "Ten patients attached, linear control"
        extends Modelica.Icons.Example;
        extends SystemModels.StandardPatientMixSystem;
        annotation(
          Diagram(coordinateSystem(extent = {{-200, -100}, {140, 100}})),
          experiment(StopTime = 15, Interval = 5e-3, Tolerance = 1e-06));
      end Scenario3;

      model Scenario4 "Nine patients attached, one more is attached at time = 10, linear control"
        extends Modelica.Icons.Example;
        extends SystemModels.StandardPatientMixSystem(
          RR(startTime = {-10, -10, -10, -10, -10, 10, -10, -10, -10, -10}));
        annotation(
          Diagram(coordinateSystem(extent = {{-200, -100}, {140, 100}})),
          experiment(StopTime = 15, Interval = 5e-3, Tolerance = 1e-06));
      end Scenario4;

      model Scenario5 "Worst-worst case, ten patients attached, all with same phase, linear control"
        extends Modelica.Icons.Example;
        extends SystemModels.WorstCasePatientSystem;
        annotation(
          Diagram(coordinateSystem(extent = {{-200, -100}, {140, 100}})),
          experiment(StopTime = 25, Interval = 5e-3, Tolerance = 1e-06));
      end Scenario5;

      model Scenario6 "Ten patients attached,leak on the last patient at time = 5"
        extends Modelica.Icons.Example;
        extends SystemModels.StandardPatientMixSystem(
          RR(each startTime = -10),
          patients(leakOpening(y = {0, 0, 0, 0, 0, 0, 0, 0, 0, if time < 5 then 0 else 1})));
        annotation(
          Diagram(coordinateSystem(extent = {{-200, -100}, {140, 100}})),
          experiment(StopTime = 15, Interval = 5e-3, Tolerance = 1e-06));
      end Scenario6;

      model Scenario7 "Twenty patients attached, linear control"
        extends Modelica.Icons.Example;
        extends SystemModels.StandardPatientMix2XSystem(
          RR(each startTime = -10));
        annotation(
          Diagram(coordinateSystem(extent = {{-200, -100}, {140, 100}})),
          experiment(StopTime = 25, Interval = 5e-3, Tolerance = 1e-06));
      end Scenario7;
    end LinearControl;

    package OnOffControl "Scenarios with on-off emergency control of the bell jar"
      extends Modelica.Icons.ExamplesPackage;

      model Scenario1 "System turned on, first patient attached at time = 10, on-off control"
        extends LinearControl.Scenario1(ambient(useOnOffControl = true));
        annotation(
          experiment(StopTime = 25, Interval = 0.005, Tolerance = 1e-06));
      end Scenario1;

      model Scenario2 "Five patients attached, on-off control"
        extends LinearControl.Scenario2(ambient(useOnOffControl = true));
        annotation(
          experiment(StopTime = 40, Interval = 0.005, Tolerance = 1e-06));
      end Scenario2;

      model Scenario3 "Ten patients attached, on-off control"
        extends LinearControl.Scenario3(ambient(useOnOffControl = true));
        annotation(
          experiment(StopTime = 25, Interval = 0.005, Tolerance = 1e-06));
      end Scenario3;

      model Scenario4 "Nine patients attached, one more is attached at time = 10, on-off control"
        extends LinearControl.Scenario4(ambient(useOnOffControl = true));
        annotation(
          experiment(StopTime = 25, Interval = 0.005, Tolerance = 1e-06));
      end Scenario4;

      model Scenario5 "Worst-worst case, ten patients attached, all with same phase, on-off control"
        extends LinearControl.Scenario5(ambient(useOnOffControl = true));
        annotation(
          experiment(StopTime = 25, Interval = 0.005, Tolerance = 1e-06));
      end Scenario5;

      model Scenario6 "Ten patients attached,leak on the last patient at time = 5"
        extends LinearControl.Scenario6(ambient(useOnOffControl = true));
        annotation(
          experiment(StopTime = 25, Interval = 0.005, Tolerance = 1e-06));
      end Scenario6;

      model Scenario7 "Twenty patients attached, on-off control"
        extends LinearControl.Scenario7(ambient(useOnOffControl = true));
        annotation(
          experiment(StopTime = 30, Interval = 0.005, Tolerance = 1e-06));
      end Scenario7;
    end OnOffControl;
  end Simulations;

  package Test "Test cases for individual components"
    extends Modelica.Icons.ExamplesPackage;

    model TestCylinder
      extends Modelica.Icons.Example;
      Components.Cylinder cylinder(A = 0.1, Vref = 0.050, yref = 0.50, ystart = 0.50, ymax = 0.80, pstart = 2400, init_y = true) annotation(
        Placement(visible = true, transformation(origin = {0, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Sources.FlowSource inletFlow annotation(
        Placement(visible = true, transformation(origin = {-30, -18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      MEV.Sources.FlowSource outletFlow annotation(
        Placement(visible = true, transformation(origin = {34, -18}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Step step1(height = 600, offset = 0, startTime = 2) annotation(
        Placement(visible = true, transformation(origin = {-72, -18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Step step2(height = -600, offset = 0, startTime = 4) annotation(
        Placement(visible = true, transformation(origin = {44, 18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      inner Ambient ambient annotation(
        Placement(transformation(extent = {{-60, 20}, {-40, 40}})));
      Modelica.Mechanics.Translational.Sources.ConstantForce constantForce(f_constant = -2400 * 3.14 * 0.2 ^ 2 / 4) annotation(
        Placement(transformation(extent = {{-28, 10}, {-8, 30}})));
    equation
      connect(inletFlow.q, step1.y) annotation(
        Line(points = {{-40, -18}, {-61, -18}}, color = {0, 0, 127}));
      connect(step2.y, outletFlow.q) annotation(
        Line(points = {{55, 18}, {74, 18}, {74, -18}, {44, -18}}, color = {0, 0, 127}));
      connect(constantForce.flange, cylinder.pistonFlange) annotation(
        Line(points = {{-8, 20}, {0, 20}, {0, -8}}, color = {0, 127, 0}));
      connect(inletFlow.port, cylinder.inlet) annotation(
        Line(points = {{-20, -18}, {-6, -18}}, color = {111, 164, 171}));
      connect(cylinder.outlet, outletFlow.port) annotation(
        Line(points = {{6, -18}, {24, -18}}, color = {111, 164, 171}));
      annotation(
        experiment(StopTime = 5, Interval = 0.06, Tolerance = 1e-06),
        Diagram(coordinateSystem(extent = {{-100, -60}, {100, 60}})),
        Icon(coordinateSystem(extent = {{-100, -60}, {100, 60}})),
        Documentation(info = "<html>
<p>Simualtion with a 70-litre tank @ 1 meter height. Initially, a flow of 10 l/s is blown into the cylinder for 2 seconds, totalling 20 liters and increasing the level from 0.5 to 0.7, then an equal outgoing flow is established, bringing the level at equilibrium.</p>
</html>"));
    end TestCylinder;

    model TestPipeSegment
      extends Modelica.Icons.Example;
      MEV.Sources.PressureSource sink annotation(
        Placement(visible = true, transformation(origin = {80, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Sources.FlowSource source annotation(
        Placement(visible = true, transformation(origin = {-30, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Step flowRate(height = 100, startTime = 0.001) annotation(
        Placement(visible = true, transformation(origin = {-70, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      ParameterizedComponents.PipeSegment pipeSegment annotation(
        Placement(transformation(extent = {{14, -10}, {34, 10}})));
      inner Ambient ambient annotation(
        Placement(transformation(extent = {{50, 30}, {70, 50}})));
    equation
      connect(flowRate.y, source.q) annotation(
        Line(points = {{-59, 0}, {-40, 0}}, color = {0, 0, 127}));
      connect(source.port, pipeSegment.inlet) annotation(
        Line(points = {{-20, 0}, {14, 0}}, color = {111, 164, 171}));
      connect(sink.port, pipeSegment.outlet) annotation(
        Line(points = {{80, 0}, {34, 0}}, color = {111, 164, 171}));
      annotation(
        experiment(StopTime = 0.002, Tolerance = 1e-06),
        Diagram(coordinateSystem(extent = {{-100, -60}, {100, 60}})),
        Icon(coordinateSystem(extent = {{-100, -60}, {100, 60}})));
    end TestPipeSegment;

    model TestPipes
      extends Modelica.Icons.Example;
      ParameterizedComponents.PipeSegment pipeSegment1(L = 2.2 * 5) annotation(
        Placement(visible = true, transformation(origin = {28, 0}, extent = {{20, -20}, {-20, 20}}, rotation = 0)));
      MEV.Sources.PressureSource sink(prel = 2400) annotation(
        Placement(visible = true, transformation(origin = {80, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Sources.FlowSource source annotation(
        Placement(visible = true, transformation(origin = {-60, -22}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
      Modelica.Blocks.Sources.Step flowRate(height = -50, startTime = 0.1) annotation(
        Placement(visible = true, transformation(origin = {-104, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Sources.FlowSource source1 annotation(
        Placement(visible = true, transformation(origin = {0, -22}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
      ParameterizedComponents.PipeSegment pipeSegment2(L = 2.2 * 5) annotation(
        Placement(visible = true, transformation(origin = {-30, 0}, extent = {{20, -20}, {-20, 20}}, rotation = 0)));
      inner Ambient ambient annotation(
        Placement(transformation(extent = {{50, 30}, {70, 50}})));
    equation
      connect(flowRate.y, source.q) annotation(
        Line(points = {{-93, -40}, {-60, -40}, {-60, -32}}, color = {0, 0, 127}));
      connect(flowRate.y, source1.q) annotation(
        Line(points = {{-93, -40}, {0, -40}, {0, -32}}, color = {0, 0, 127}));
      connect(pipeSegment1.inlet, sink.port) annotation(
        Line(points = {{48, 0}, {80, 0}}, color = {111, 164, 171}));
      connect(pipeSegment2.inlet, pipeSegment1.outlet) annotation(
        Line(points = {{-10, 0}, {8, 0}}, color = {111, 164, 171}));
      connect(source.port, pipeSegment2.outlet) annotation(
        Line(points = {{-60, -12}, {-60, 0}, {-50, 0}}, color = {111, 164, 171}));
      connect(source1.port, pipeSegment1.outlet) annotation(
        Line(points = {{0, -12}, {0, 0}, {8, 0}}, color = {111, 164, 171}));
      annotation(
        experiment(StopTime = 10, Tolerance = 1e-06),
        Diagram(coordinateSystem(extent = {{-120, -60}, {120, 60}})),
        Icon(coordinateSystem(extent = {{-120, -60}, {120, 60}})));
    end TestPipes;

    model TestPatientValve
      extends Modelica.Icons.Example;
      ParameterizedComponents.PatientValve valve annotation(
        Placement(visible = true, transformation(origin = {0, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      MEV.Sources.PressureSource sink annotation(
        Placement(visible = true, transformation(origin = {38, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      MEV.Sources.PressureSource source(prel = 1200) annotation(
        Placement(visible = true, transformation(origin = {-38, -10}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Step opening(startTime = 1) annotation(
        Placement(visible = true, transformation(origin = {-22, 22}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      inner Ambient ambient annotation(
        Placement(transformation(extent = {{40, 40}, {60, 60}})));
    equation
      connect(opening.y, valve.opening) annotation(
        Line(points = {{-11, 22}, {0, 22}, {0, -5.8}, {0, -5.8}}, color = {0, 0, 127}));
      connect(valve.outlet, sink.port) annotation(
        Line(points = {{10, -10}, {38, -10}}, color = {0, 0, 255}));
      connect(source.port, valve.inlet) annotation(
        Line(points = {{-38, -10}, {-10, -10}, {-10, -10}, {-10, -10}}, color = {0, 0, 255}));
      annotation(
        experiment(StartTime = 0, StopTime = 2, Tolerance = 1e-6, Interval = 0.004));
    end TestPatientValve;

    model TestSupplyValve
      extends Modelica.Icons.Example;
      ParameterizedComponents.SupplyValve valve annotation(
        Placement(visible = true, transformation(origin = {0, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      MEV.Sources.PressureSource sink annotation(
        Placement(visible = true, transformation(origin = {38, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      MEV.Sources.PressureSource source(prel = 2e5, useRelativePressure = true) annotation(
        Placement(visible = true, transformation(origin = {-38, -10}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Step opening(startTime = 1) annotation(
        Placement(visible = true, transformation(origin = {-22, 22}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      inner Ambient ambient annotation(
        Placement(transformation(extent = {{40, 40}, {60, 60}})));
    equation
      connect(opening.y, valve.opening) annotation(
        Line(points = {{-11, 22}, {0, 22}, {0, -5.8}, {0, -5.8}}, color = {0, 0, 127}));
      connect(valve.outlet, sink.port) annotation(
        Line(points = {{10, -10}, {38, -10}}, color = {0, 0, 255}));
      connect(source.port, valve.inlet) annotation(
        Line(points = {{-38, -10}, {-10, -10}, {-10, -10}, {-10, -10}}, color = {0, 0, 255}));
      annotation(
        experiment(StartTime = 0, StopTime = 2, Tolerance = 1e-6, Interval = 0.004));
    end TestSupplyValve;

    model TestStandardPatient
      extends Modelica.Icons.Example;
      MEV.Sources.PressureSource source(prel = 2400) annotation(
        Placement(visible = true, transformation(origin = {0, 34}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Pulse opening(width = 40, period = 60 / 20, startTime = 0) annotation(
        Placement(visible = true, transformation(origin = {-54, -16}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      ParameterizedComponents.StandardPatient standardPatient annotation(
        Placement(transformation(extent = {{-20, -36}, {20, 4}})));
      inner Ambient ambient annotation(
        Placement(transformation(extent = {{40, 60}, {60, 80}})));
    equation
      connect(source.port, standardPatient.supply) annotation(
        Line(points = {{0, 34}, {0, 0}}, color = {111, 164, 171}));
      connect(standardPatient.valveCommand, opening.y) annotation(
        Line(points = {{-16.8, -16}, {-43, -16}}, color = {0, 0, 127}));
      annotation(
        experiment(StopTime = 10, Tolerance = 1e-06),
        Diagram(coordinateSystem(extent = {{-120, -120}, {120, 120}})),
        Icon(coordinateSystem(extent = {{-120, -120}, {120, 120}})));
    end TestStandardPatient;

    model TestBellJarOpenLoop
      extends Modelica.Icons.Example;
      ParameterizedComponents.BellJar bellJar(A = 3.14 * 10e-2 ^ 2, ystart = 0.8, ymax = 1.2, pdes = 2400, d = 10, pstart = 2400) annotation(
        Placement(transformation(extent = {{-16, -24}, {16, 16}})));
      Sources.PressureSource pressureSource(prel = 2e5, useRelativePressure = true) annotation(
        Placement(transformation(extent = {{-130, -30}, {-110, -10}})));
      ParameterizedComponents.SupplyValve supplyValve(Av = 1.4e-5, wnom = 0.001, dpnom = 2e5) annotation(
        Placement(transformation(extent = {{-68, -28}, {-52, -12}})));
      inner Ambient ambient annotation(
        Placement(transformation(extent = {{110, 30}, {130, 50}})));
      Modelica.Blocks.Sources.Step step(height = 0.1, startTime = 1) annotation(
        Placement(transformation(extent = {{-100, 0}, {-86, 14}})));
    equation
      connect(pressureSource.port, supplyValve.inlet) annotation(
        Line(points = {{-120, -20}, {-68, -20}}, color = {111, 164, 171}));
      connect(supplyValve.outlet, bellJar.inlet) annotation(
        Line(points = {{-52, -20}, {-12, -20}}, color = {111, 164, 171}));
      connect(step.y, supplyValve.opening) annotation(
        Line(points = {{-85.3, 7}, {-60, 7}, {-60, -16.64}}, color = {0, 0, 127}));
      annotation(
        Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-140, -60}, {140, 60}})),
        Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-140, -60}, {140, 60}})),
        experiment(StopTime = 10, Interval = 0.005, Tolerance = 1e-06));
    end TestBellJarOpenLoop;

    model TestBellJarLinearControl
      extends Modelica.Icons.Example;
      Components.ControlValve controlValve(Av = 2.5e-6, wnom = 0.001, dpnom = 2e5) annotation(
        Placement(transformation(extent = {{-68, -28}, {-52, -12}})));
      Sources.PressureSource pressureSource(prel = 3e5, useRelativePressure = true) annotation(
        Placement(transformation(extent = {{-130, -30}, {-110, -10}})));
      Sources.FlowSource flowSource annotation(
        Placement(transformation(extent = {{60, -30}, {40, -10}})));
      Modelica.Blocks.Sources.Pulse opening(amplitude = -25, width = 40, period = 60 / 20, startTime = 0) annotation(
        Placement(visible = true, transformation(origin = {52, 12}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      inner Ambient ambient annotation(
        Placement(transformation(extent = {{110, 30}, {130, 50}})));
      Controls.LinearController controller(level_sp = (bellJar.ymax + bellJar.ymin) / 2, Kp = 10, T = 0.3) annotation(
        Placement(transformation(extent = {{-26, 0}, {-42, 16}})));
      ParameterizedComponents.BellJar bellJar(ystart = 0.06) annotation(
        Placement(transformation(extent = {{-12, -24}, {20, 16}})));
    equation
      connect(pressureSource.port, controlValve.inlet) annotation(
        Line(points = {{-120, -20}, {-68, -20}}, color = {111, 164, 171}));
      connect(opening.y, flowSource.q) annotation(
        Line(points = {{63, 12}, {80, 12}, {80, -20}, {60, -20}}, color = {0, 0, 127}));
      connect(controller.valveOpening, controlValve.opening) annotation(
        Line(points = {{-42.8, 8}, {-60, 8}, {-60, -16.64}}, color = {0, 0, 127}));
      connect(bellJar.inlet, controlValve.outlet) annotation(
        Line(points = {{-8, -20}, {-52, -20}}, color = {111, 164, 171}));
      connect(bellJar.outlet, flowSource.port) annotation(
        Line(points = {{16, -20}, {40, -20}}, color = {111, 164, 171}));
      connect(bellJar.y, controller.level) annotation(
        Line(points = {{-12.4, 8}, {-24.4, 8}}, color = {0, 0, 127}));
      annotation(
        Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-140, -60}, {140, 60}})),
        Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-140, -60}, {140, 60}})),
        experiment(StopTime = 100, StartTime = 0, Tolerance = 1e-06, Interval = 0.2));
    end TestBellJarLinearControl;

    model TestCylinderUnitOnOff
      extends Modelica.Icons.Example;
      Sources.PressureSource pressureSource(prel = 3e5, useRelativePressure = true) annotation(
        Placement(transformation(extent = {{-130, -30}, {-110, -10}})));
      Components.ControlValve controlValve(Av = 2.5e-6, wnom = 0.001, dpnom = 2e5) annotation(
        Placement(transformation(extent = {{-68, -28}, {-52, -12}})));
      Sources.FlowSource flowSource annotation(
        Placement(transformation(extent = {{60, -30}, {40, -10}})));
      Modelica.Blocks.Sources.Pulse opening(amplitude = -25, width = 40, period = 60 / 20, startTime = 0) annotation(
        Placement(visible = true, transformation(origin = {52, 12}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      inner Ambient ambient(useOnOffControl = true) annotation(
        Placement(transformation(extent = {{110, 30}, {130, 50}})));
      Controls.OnOffControllerWithHysteresis controller(level_min = bellJar.ymin + 0.02, level_max = bellJar.ymax - 0.02) annotation(
        Placement(transformation(extent = {{-28, 0}, {-44, 16}})));
      ParameterizedComponents.BellJar bellJar annotation(
        Placement(transformation(extent = {{-10, -24}, {22, 16}})));
    equation
      connect(pressureSource.port, controlValve.inlet) annotation(
        Line(points = {{-120, -20}, {-68, -20}}, color = {111, 164, 171}));
      connect(opening.y, flowSource.q) annotation(
        Line(points = {{63, 12}, {80, 12}, {80, -20}, {60, -20}}, color = {0, 0, 127}));
      connect(controller.valveOpening, controlValve.opening) annotation(
        Line(points = {{-44.8, 8}, {-60, 8}, {-60, -16.64}}, color = {0, 0, 127}));
      connect(flowSource.port, bellJar.outlet) annotation(
        Line(points = {{40, -20}, {18, -20}}, color = {111, 164, 171}));
      connect(controlValve.outlet, bellJar.inlet) annotation(
        Line(points = {{-52, -20}, {-6, -20}}, color = {111, 164, 171}));
      connect(bellJar.y, controller.level) annotation(
        Line(points = {{-10.4, 8}, {-26.4, 8}}, color = {0, 0, 127}));
      annotation(
        Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-140, -60}, {140, 60}})),
        Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-140, -60}, {140, 60}})),
        experiment(StopTime = 100, StartTime = 0, Tolerance = 1e-06, Interval = 0.2));
    end TestCylinderUnitOnOff;

    model TestCylinderUnitLinearize
      extends Modelica.Icons.Example;
      Sources.PressureSource pressureSource(prel = 2e5, useRelativePressure = true) annotation(
        Placement(transformation(extent = {{-130, -30}, {-110, -10}})));
      ParameterizedComponents.SupplyValve supplyValve annotation(
        Placement(transformation(extent = {{-68, -28}, {-52, -12}})));
      inner Ambient ambient annotation(
        Placement(transformation(extent = {{110, 30}, {130, 50}})));
      Modelica.Blocks.Interfaces.RealInput opening annotation(
        Placement(transformation(extent = {{-160, 8}, {-120, 48}}), iconTransformation(extent = {{-120, -20}, {-80, 20}})));
      Modelica.Blocks.Interfaces.RealOutput y "Height of piston [m]" annotation(
        Placement(transformation(extent = {{-30, 18}, {-10, 38}}), iconTransformation(extent = {{90, -10}, {110, 10}})));
      ParameterizedComponents.BellJar bellJar annotation(
        Placement(transformation(extent = {{-28, -24}, {4, 16}})));
    equation
      connect(pressureSource.port, supplyValve.inlet) annotation(
        Line(points = {{-120, -20}, {-68, -20}}, color = {111, 164, 171}));
      connect(supplyValve.opening, opening) annotation(
        Line(points = {{-60, -16.64}, {-60, 28}, {-140, 28}}, color = {0, 0, 127}));
      connect(supplyValve.outlet, bellJar.inlet) annotation(
        Line(points = {{-52, -20}, {-24, -20}}, color = {111, 164, 171}));
      connect(bellJar.y, y) annotation(
        Line(points = {{-28.4, 8}, {-42, 8}, {-42, 28}, {-20, 28}}, color = {0, 0, 127}));
      annotation(
        Icon(coordinateSystem(preserveAspectRatio = false)),
        Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-140, -60}, {140, 60}})),
        experiment(StopTime = 10, Interval = 0.005, Tolerance = 1e-06));
    end TestCylinderUnitLinearize;

    model TestDutyCycleGenerator
      extends Modelica.Icons.Example;
      Controls.DutyCycleGenerator dutyCycleGenerator(T = 3, phase = 180) annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Step step(height = 20, offset = 50, startTime = 10) annotation(
        Placement(visible = true, transformation(origin = {-46, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(step.y, dutyCycleGenerator.dutyCycle) annotation(
        Line(points = {{-34, 0}, {-10, 0}, {-10, 0}, {-10, 0}}, color = {0, 0, 127}));
      annotation(
        experiment(StartTime = 0, StopTime = 20, Tolerance = 1e-06, Interval = 0.04));
    end TestDutyCycleGenerator;

    model TestSequencer
      extends Modelica.Icons.Example;
      parameter Integer N = 20;
      Real sums[N];
      Controls.PatientValveSequencer sequencer(N = 20) annotation(
        Placement(transformation(extent = {{-10, -12}, {14, 12}})));
      Modelica.Blocks.Sources.Step step[N](height = fill(50, N)) annotation(
        Placement(transformation(extent = {{-58, -10}, {-38, 10}})));
    equation
      for i in 1:N loop
        sums[i] = sum(sequencer.valveCommands[j] for j in 1:i);
      end for;
      connect(step.y, sequencer.dutyCycles) annotation(
        Line(points = {{-37, 0}, {-10, 0}}, color = {0, 0, 127}));
      annotation(
        experiment(StopTime = 10, StartTime = 0, Tolerance = 1e-06, Interval = 0.02),
        Icon(coordinateSystem(preserveAspectRatio = false)),
        Diagram(coordinateSystem(preserveAspectRatio = false)));
    end TestSequencer;

    model TestPatientController1 "Test patient controller with initial RR > 1"
      extends Modelica.Icons.Example;
      Controls.PatientController patientController annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
      Modelica.Blocks.Sources.TimeTable RR(table = [0, 30; 48, 30; 48, 20; 79, 20; 79, 0; 100, 0]) annotation(
        Placement(visible = true, transformation(origin = {-62, 12}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.TimeTable dutyCycle(table = [0, 50; 20, 50; 50, 25; 100, 25]) annotation(
        Placement(visible = true, transformation(origin = {-62, -24}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(RR.y, patientController.RR) annotation(
        Line(points = {{-50, 12}, {-28, 12}, {-28, 6}, {-12, 6}}, color = {0, 0, 127}));
      connect(dutyCycle.y, patientController.dutyCycle) annotation(
        Line(points = {{-50, -24}, {-36, -24}, {-36, -6}, {-12, -6}}, color = {0, 0, 127}));
      annotation(
        experiment(StartTime = 0, StopTime = 100, Tolerance = 1e-06, Interval = 0.2));
    end TestPatientController1;

    model TestPatientController2
      extends MEV.Test.TestPatientController1(
        RR(table = [0, 0; 10, 0; 10, 30; 48, 30; 48, 20; 79, 20; 79, 0; 100, 0]));
      annotation(
        experiment(StartTime = 0, StopTime = 100, Tolerance = 1e-06, Interval = 0.2));
    end TestPatientController2;
  end Test;

  package Media "Medium models"
    extends Modelica.Icons.MaterialPropertiesPackage;

    model AirAdiabatic "Air model assuming adiabatic reversible compression"
      outer Ambient ambient "Ambient model";
      constant Modelica.SIunits.SpecificHeatCapacity Rstar = ambient.Rstar "Gas constant for air";
      constant Modelica.SIunits.MolarMass MM = ambient.MM "Molar mass of air";
      constant Modelica.SIunits.PerUnit gamma = 1.4 "cp/cv";
      connector InputAbsolutePressure = input Types.AbsolutePressure;
      Types.Temperature Tref = ambient.T "Reference gas temperature";
      Types.AbsolutePressure pref = ambient.p "Reference gas pressure";
      InputAbsolutePressure p "Absolute pressure";
      Types.Density rho "Density";
    equation
      rho = pref / (Rstar * Tref) * (1 + (p - pref) / pref / gamma);
      annotation(
        Icon(coordinateSystem(preserveAspectRatio = false)),
        Diagram(coordinateSystem(preserveAspectRatio = false)));
    end AirAdiabatic;

    model Air "Air model"
      outer Ambient ambient "Ambient model";
      constant Modelica.SIunits.SpecificHeatCapacity Rstar = ambient.Rstar "Gas constant for air";
      constant Modelica.SIunits.MolarMass MM = ambient.MM "Molar mass of air";
      connector InputAbsolutePressure = input Types.AbsolutePressure;
      Types.Temperature Tref = ambient.T "Reference gas temperature";
      InputAbsolutePressure p "Absolute pressure";
      Types.Density rho "Density";
    equation
      rho = p / (Rstar * Tref);
      annotation(
        Icon(coordinateSystem(preserveAspectRatio = false)),
        Diagram(coordinateSystem(preserveAspectRatio = false)));
    end Air;
    annotation(
      Icon(coordinateSystem(preserveAspectRatio = false)),
      Diagram(coordinateSystem(preserveAspectRatio = false)));
  end Media;

  package Types "Custom types definitions"
    extends Modelica.Icons.TypesPackage;
    type RelativePressure = Modelica.SIunits.PressureDifference(nominal = 1000, displayUnit = "Pa");
    type AbsolutePressure = Modelica.SIunits.Pressure(nominal = 1e5, displayUnit = "bar");
    type MassFlowRate = Modelica.SIunits.MassFlowRate(nominal = 1e-4);
    type VolumeFlowRate = Modelica.SIunits.VolumeFlowRate(nominal = 1e-4);
    type Volume = Modelica.SIunits.Volume(displayUnit = "l");
    type Density = Modelica.SIunits.Density(nominal = 1.2, start = 1.2, displayUnit = "kg/m3");
    type Length = Modelica.SIunits.Length;
    type Velocity = Modelica.SIunits.Velocity;
    type VelocitySquared = Real(unit = "m2/s2");
    type Temperature = Modelica.SIunits.Temperature(start = 288.15, nominal = 300);
    type Compliance = Real(unit = "m3/Pa");
    type Resistance = Real(unit = "Pa.s/m3");
  end Types;

  package Controls
    extends Modelica.Icons.Package;

    model OnOffControllerWithHysteresis
      extends Modelica.Blocks.Icons.Block;
      outer Ambient ambient "Ambient model";
      Modelica.Blocks.Interfaces.RealInput level "Level of bell jar [m]" annotation(
        Placement(transformation(extent = {{-140, -20}, {-100, 20}})));
      Modelica.Blocks.Interfaces.RealOutput valveOpening "Valve opening in p.u." annotation(
        Placement(transformation(extent = {{100, -10}, {120, 10}})));
      parameter Types.Length level_min "Minimum level threshold";
      parameter Types.Length level_max "Maximum level threshold";
      Boolean on;
    initial equation
      if ambient.useOnOffControl then
        on = if level < level_max then true else false;
      else
        on = false;
      end if;
    equation
      when level > level_max then
        on = false;
      elsewhen level < level_min and ambient.useOnOffControl then
        on = true;
      end when;
      valveOpening = if on then 1 else 0;
      annotation(
        Icon(coordinateSystem(preserveAspectRatio = false), graphics = {Text(extent = {{-78, 52}, {74, -44}}, lineColor = {28, 108, 200}, textString = "On-Off")}),
        Diagram(coordinateSystem(preserveAspectRatio = false)));
    end OnOffControllerWithHysteresis;

    model LinearController
      extends Modelica.Blocks.Icons.Block;
      outer Ambient ambient "Ambient model";
      parameter Types.Length level_sp "Level set point";
      parameter Real Kp "Proportional Gain";
      parameter Modelica.SIunits.Time T "Time Constant";
      Modelica.Blocks.Interfaces.RealInput level "Level of bell jar [m]" annotation(
        Placement(transformation(extent = {{-20, -20}, {20, 20}}, rotation = 90, origin = {-36, -64}), iconTransformation(extent = {{-140, -20}, {-100, 20}})));
      Modelica.Blocks.Interfaces.RealOutput valveOpening "Valve opening in p.u." annotation(
        Placement(transformation(extent = {{100, -10}, {120, 10}})));
      Modelica.Blocks.Continuous.FirstOrder firstOrder(k = Kp, T = T, initType = Modelica.Blocks.Types.Init.SteadyState) annotation(
        Placement(transformation(extent = {{-8, -18}, {12, 2}})));
      Modelica.Blocks.Math.Feedback feedback annotation(
        Placement(transformation(extent = {{-46, -18}, {-26, 2}})));
      Modelica.Blocks.Sources.RealExpression setPoint(y = level_sp) annotation(
        Placement(transformation(extent = {{-82, -18}, {-62, 2}})));
      Modelica.Blocks.Logical.Switch switch1 annotation(
        Placement(transformation(extent = {{58, -10}, {78, 10}})));
      Modelica.Blocks.Sources.RealExpression zero annotation(
        Placement(transformation(extent = {{-2, 38}, {24, 66}})));
      Modelica.Blocks.Sources.BooleanExpression useOnOffControl(y = ambient.useOnOffControl) annotation(
        Placement(transformation(extent = {{-36, 14}, {26, 40}})));
    equation
      connect(feedback.u2, level) annotation(
        Line(points = {{-36, -16}, {-36, -64}}, color = {0, 0, 127}));
      connect(feedback.y, firstOrder.u) annotation(
        Line(points = {{-27, -8}, {-10, -8}}, color = {0, 0, 127}));
      connect(setPoint.y, feedback.u1) annotation(
        Line(points = {{-61, -8}, {-44, -8}}, color = {0, 0, 127}));
      connect(firstOrder.y, switch1.u3) annotation(
        Line(points = {{13, -8}, {56, -8}}, color = {0, 0, 127}));
      connect(switch1.y, valveOpening) annotation(
        Line(points = {{79, 0}, {110, 0}}, color = {0, 0, 127}));
      connect(zero.y, switch1.u1) annotation(
        Line(points = {{25.3, 52}, {44, 52}, {44, 8}, {56, 8}}, color = {0, 0, 127}));
      connect(switch1.u2, useOnOffControl.y) annotation(
        Line(points = {{56, 0}, {36, 0}, {36, 27}, {29.1, 27}}, color = {255, 0, 255}));
      annotation(
        Icon(coordinateSystem(preserveAspectRatio = false), graphics = {Text(extent = {{-78, 52}, {74, -44}}, lineColor = {28, 108, 200}, textString = "Linear")}),
        Diagram(coordinateSystem(preserveAspectRatio = false)));
    end LinearController;

    model PatientValveSequencer
      parameter Integer N = 10 "Max number of patients";
      parameter Real Nr = 20 "Number of respiratory acts per minute";
      parameter Real phases[:] = {0, 180, 90, 270, 45, 225, 135, 315, 22, 202, 112, 292, 158, 338, 68, 158, 11, 191, 281, 56, 237, 146} "Phase angle of each patient";
      final parameter Modelica.SIunits.Time T = 60 / Nr;
      Modelica.Blocks.Interfaces.RealOutput valveCommands[N] annotation(
        Placement(transformation(extent = {{-338, 6}, {-318, 26}}), iconTransformation(extent = {{100, -20}, {140, 20}})));
      Modelica.Blocks.Interfaces.RealInput dutyCycles[N] "Duty cycle of each patient in percent" annotation(
        Placement(transformation(extent = {{-282, -4}, {-242, 36}}), iconTransformation(extent = {{-140, -20}, {-100, 20}})));
      MEV.Controls.DutyCycleGenerator generators[N](each T = T, phase = phases[1:N]);
    equation
      connect(generators.dutyCycle, dutyCycles);
      connect(generators.opening, valveCommands);
      annotation(
        Icon(coordinateSystem(extent = {{-120, -120}, {120, 120}}), graphics = {Rectangle(extent = {{-120, 120}, {120, -120}}, lineColor = {0, 0, 0}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid), Ellipse(extent = {{-80, 80}, {80, -80}}, lineColor = {160, 160, 164}), Line(points = {{0, 80}, {0, 60}}, color = {160, 160, 164}), Line(points = {{80, 0}, {60, 0}}, color = {160, 160, 164}), Line(points = {{0, -80}, {0, -60}}, color = {160, 160, 164}), Line(points = {{-80, 0}, {-60, 0}}, color = {160, 160, 164}), Line(points = {{37, 70}, {26, 50}}, color = {160, 160, 164}), Line(points = {{70, 38}, {49, 26}}, color = {160, 160, 164}), Line(points = {{71, -37}, {52, -27}}, color = {160, 160, 164}), Line(points = {{39, -70}, {29, -51}}, color = {160, 160, 164}), Line(points = {{-39, -70}, {-29, -52}}, color = {160, 160, 164}), Line(points = {{-71, -37}, {-50, -26}}, color = {160, 160, 164}), Line(points = {{-71, 37}, {-54, 28}}, color = {160, 160, 164}), Line(points = {{-38, 70}, {-28, 51}}, color = {160, 160, 164}), Line(points = {{0, 0}, {44, 48}}, thickness = 0.5)}),
        Diagram(coordinateSystem(extent = {{-120, -120}, {120, 120}})));
    end PatientValveSequencer;

    model DutyCycleGenerator
      parameter Modelica.SIunits.Time T "Sampling time";
      parameter Real phase "Phase angle in degrees";
      Modelica.Blocks.Interfaces.RealInput dutyCycle annotation(
        Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-100, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealOutput opening annotation(
        Placement(visible = true, transformation(origin = {102, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      discrete Modelica.SIunits.Time offTime(start = -1e6, fixed = true);
    equation
      when sample(phase / 360 * T, T) then
        offTime = time + T * dutyCycle / 100;
      end when;
      opening = if time < offTime then 1 else 0;
      annotation(
        Icon(graphics = {Rectangle(extent = {{-100, 100}, {100, -100}})}));
    end DutyCycleGenerator;

    model PatientController
      Modelica.Blocks.Interfaces.RealOutput opening "Opening signal for patient valve [p.u.]" annotation(
        Placement(visible = true, transformation(extent = {{84, -2}, {104, 18}}, rotation = 0), iconTransformation(extent = {{80, -20}, {120, 20}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput dutyCycle "Duty cycle [%]" annotation(
        Placement(visible = true, transformation(extent = {{-130, -36}, {-90, 4}}, rotation = 0), iconTransformation(extent = {{-120, -80}, {-80, -40}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput RR "Respiratory rate [acts/min]" annotation(
        Placement(visible = true, transformation(extent = {{-130, 16}, {-90, 56}}, rotation = 0), iconTransformation(extent = {{-120, 40}, {-80, 80}}, rotation = 0)));
      discrete Modelica.SIunits.Time t_on "Next opening time";
      discrete Modelica.SIunits.Time t_off "Next closing time";
      Boolean active "The controller is activated";
    initial equation
      active = if RR > 1 then true else false;
      if active then
        t_on = time + 60/RR;
        t_off = time + 60 / RR * dutyCycle / 100;
      else
        t_on = time;
        t_off = time;
      end if;
    equation
      when RR > 1 then
        active = true;
      elsewhen RR < 1 then
        active = false;
      end when;
      when pre(active) and time >= pre(t_on) then
        t_on = time + 60 / RR;
        t_off = time + 60 / RR * dutyCycle / 100;
      end when;
      opening = if active and time < t_off then 1 else 0;
      annotation(
        Icon(coordinateSystem(extent = {{-120, -120}, {120, 120}}), graphics = {Rectangle(fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}), Ellipse(lineColor = {160, 160, 164}, extent = {{-80, 80}, {80, -80}}, endAngle = 360), Line(points = {{0, 80}, {0, 60}}, color = {160, 160, 164}), Line(points = {{80, 0}, {60, 0}}, color = {160, 160, 164}), Line(points = {{0, -80}, {0, -60}}, color = {160, 160, 164}), Line(points = {{-80, 0}, {-60, 0}}, color = {160, 160, 164}), Line(points = {{37, 70}, {26, 50}}, color = {160, 160, 164}), Line(points = {{70, 38}, {49, 26}}, color = {160, 160, 164}), Line(points = {{71, -37}, {52, -27}}, color = {160, 160, 164}), Line(points = {{39, -70}, {29, -51}}, color = {160, 160, 164}), Line(points = {{-39, -70}, {-29, -52}}, color = {160, 160, 164}), Line(points = {{-71, -37}, {-50, -26}}, color = {160, 160, 164}), Line(points = {{-71, 37}, {-54, 28}}, color = {160, 160, 164}), Line(points = {{-38, 70}, {-28, 51}}, color = {160, 160, 164}), Line(points = {{0, 0}, {44, 48}}, thickness = 0.5)}),
        Diagram(coordinateSystem(extent = {{-120, -120}, {120, 120}})),
        experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-6, Interval = 0.002));
    end PatientController;
  end Controls;
  annotation(
    version = "1.1.0",
    uses(Modelica(version = "3.2.3")));
end MEV;
