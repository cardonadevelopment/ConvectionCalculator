/*******************************************************************************
 * The MIT License (MIT)
 *
 * Convection calculator app
 * Copyright (C) 2017 Luis F. Cardona S.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 ******************************************************************************/

package co.edu.itm.www.convectioncalculator;

import android.content.Intent;
import android.content.SharedPreferences;
import android.support.v7.app.AlertDialog;
import android.support.v7.preference.PreferenceManager;
import android.os.Bundle;
import android.support.v7.widget.Toolbar;
import android.view.LayoutInflater;
import android.view.View;
import android.webkit.WebView;
import android.widget.AdapterView;
import android.widget.ArrayAdapter;
import android.widget.ImageView;
import android.widget.SeekBar;
import android.widget.TextView;
import android.support.v7.app.AppCompatActivity;
import android.widget.Spinner;
import android.view.Menu;
import android.view.MenuItem;
import android.widget.Toast;

import java.text.NumberFormat;


public class MainActivity extends AppCompatActivity {

    // Used to load the 'native-lib' library on application startup.
    static {
        System.loadLibrary("CoolProp");
    }

    public String tubeGeometry = "Annular";
    public String tubeMaterial = "Steel";
    public String selectedFluid = "Water";
    public String flowInputKnownData = "massFlow";
    public String thermalInputKnownData = "surfaceTemperature";
    public String flowRegime = "Laminar";

    private TextView sALabelTextview;
    private TextView sBLabelTextview;
    private TextView sATextview;
    private TextView sBTextview;
    private TextView tempinletTextView;
    private TextView flowInputDataTextView;
    private TextView thermalInputDataTextView;
    private TextView tubeLTextView;
    private TextView diamLabelTextview;
    private TextView innerDiamLabelTextview;

    private TextView massFlowLabelsTextView;
    private TextView volumetricFlowLabelsTextView;
    private TextView inletVelocityLabelsTextView;



    //---------Results data text views------------
    private TextView outletTempDataTextView;
    private TextView massFlowDataTextView;
    private TextView volumetricFlowDataTextView;
    private TextView inletVelocityDataTextView;
    private TextView ReynoldsDataTextView;
    private TextView flowRegimeDataTextView;
    private TextView frictionFactorDataTextView;
    private TextView thermalEntryLengthDataTextView;
    private TextView hydrodynamicEntryLengthDataTextView;
    private TextView NusseltDataTextView;
    private TextView convectionCoefficientDataTextView;
    private TextView heatTransferRateDataTextView;
    private TextView NTUDataTextView;
    private TextView meanTempDataTextView;
    private TextView fluidDensityDataTextView;
    private TextView fluidCpDataTextView;
    private TextView fluidPrandtlDataTextView;
    private TextView fluidConductivityDataTextView;
    private TextView fluidKinematicViscosityDataTextView;
    private TextView fluidDynamicViscosityDataTextView;
    private TextView fluidDynamicViscositySurfaceTempDataTextView;



    //------- Fluid Properties-------------
    public double fluidDensity = 1;//kg/m3
    public double fluidConductivity = 1;//W/(m°C)
    public double fluidDynamicViscosity = 1;//kg/(m-s)
    public double fluidDynamicViscositySurfaceTemp = 1;//kg/(m-s)
    public double fluidKinematicViscosity = 1;//m2/s
    public double fluidPrandtl = 1;//non-dimensional
    public double fluidCp = 1;//J/(kg-K)
    public double inletPressure = 101300;//Pa

    public double surfaceArea = 1;//m2
    public double tubeLength= 1;//m
    public double Diameter = 2;//m
    public double innerDiameter = 0.000001;//m
    public double SideA = 1;//m
    public double SideB = 1;//m
    public double ratioAB = 1;//non-dimensional
    public double ratioDiameters = 1;
    public double sectionArea = 1;//m2
    public double inletVelocity = 1;//m/s
    public double massFlow = 1;//kg/s
    public double volumetricFlow = 1;//m3/s
    public double Reynolds = 1;//non-dimensional
    public double Nusselt = 1;//non-dimensional
    public double convectionCoefficient = 1;//W/(m2°C)
    public double heatTransferRate =1;//W
    public double NTU =1;//non-dimensional
    public double deltaTlm = 10;//°C
    public double frictionFactor = 1;//non-dimensional
    public double roughness = 0;//mm
    public double thermalEntryLength=0;//m
    public double hydrodynamicEntryLength=0;//m
    public double auxvariable = 1;
    public double auxvariable2 = 1;
    public double auxvariable3 = 1;

    private double inletTemp = 20;//°C
    private double outletTemp = 21;//°C
    private double meanTemp = 20;//°C
    private double surfaceTemp = 20;//°C


    //----Limits of variables involved in seekbars-------
    public double minInletTemperature = 1;//°C
    public double maxInletTemperature = 100;//°C
    public double minPipeLength = 0.01;//m
    public double maxPipeLength = 100;//m
    public double minSide = 0.01;//m
    public double maxSide = 1;//m
    public double minMassFlow = 0.01;//kg/s
    public double maxMassFlow = 1;//kg/s
    public double minVolumetricFlow = 0.01;//m3/s
    public double maxVolumetricFlow = 1;//m3/s
    public double minInletVelocity = 0.01;//m/s
    public double maxInletVelocity = 1;//m/s
    public double minSurfaceTemperature = 0;//°C
    public double maxSurfaceTemperature = 10;//°C
    public double minHeatTransferRate=0;//W
    public double maxHeatTransferRate=10;//W


    private double hydraulicDiameter = 1;///m
    private double flowInputKnownDataValue = 1;//this can be °C or m2 depending on the case
    private double thermalInputKnownDataValue = 1;//this can be °C or W depending on the case

    private static final NumberFormat tempFormat = NumberFormat.getInstance();
    private static final NumberFormat lengthFormat = NumberFormat.getInstance();
    private static final NumberFormat frictionFormat = NumberFormat.getInstance();
    private static final NumberFormat viscosityFormat = NumberFormat.getInstance();

    public void checkDiameters(){

        Spinner spinnerDiam = (Spinner) findViewById(R.id.DiameterSpinner);
        if (tubeGeometry == "Annular")
        {
            if(Diameter<=innerDiameter)
            {
                spinnerDiam.setSelection(15);//15 is the largest
                Diameter= 7.967*25.4/1000;//conversion to m
                Toast.makeText(this, R.string.errorMessageDiameters, Toast.LENGTH_SHORT).show();
            }
        }
    }


    public void setPropertiesVisibility(){

        meanTempDataTextView = (TextView) findViewById(R.id.meanTempTextView);
        fluidDensityDataTextView = (TextView) findViewById(R.id.fluidDensityTextView);
        fluidCpDataTextView = (TextView) findViewById(R.id.fluidCpTextView);
        fluidPrandtlDataTextView = (TextView) findViewById(R.id.fluidPrandtlTextView);
        fluidConductivityDataTextView = (TextView) findViewById(R.id.fluidConductivityTextView);
        fluidKinematicViscosityDataTextView = (TextView) findViewById(R.id.fluidKinematicViscosityTextView);
        fluidDynamicViscosityDataTextView = (TextView) findViewById(R.id.fluidDynamicViscosityTextView);
        fluidDynamicViscositySurfaceTempDataTextView = (TextView) findViewById(R.id.fluidDynamicViscositySurfaceTempTextView);
        TextView meanTempLabelTextView = (TextView) findViewById(R.id.meanTempLabelTextView);
        TextView fluidDensityLabelTextView = (TextView) findViewById(R.id.fluidDensityLabelTextView);
        TextView fluidCpLabelTextView = (TextView) findViewById(R.id.fluidCpLabelTextView);
        TextView fluidPrandtlLabelTextView = (TextView) findViewById(R.id.fluidPrandtlLabelTextView);
        TextView fluidConductivityLabelTextView = (TextView) findViewById(R.id.fluidConductivityLabelTextView);
        TextView fluidKinematicViscosityLabelTextView = (TextView) findViewById(R.id.fluidKinematicViscosityLabelTextView);
        TextView fluidDynamicViscosityLabelTextView = (TextView) findViewById(R.id.fluidDynamicViscosityLabelTextView);
        TextView fluidDynamicViscositySurfaceTempLabelTextView = (TextView) findViewById(R.id.fluidDynamicViscositySurfaceTempLabelTextView);


        SharedPreferences sharedPref = PreferenceManager.getDefaultSharedPreferences(this);
        double propertiesVisualization = Double.parseDouble(sharedPref.getString("visualizationList", "0"));

        if(propertiesVisualization==0)
        {
            meanTempDataTextView.setVisibility(View.GONE);
            fluidDensityDataTextView.setVisibility(View.GONE);
            fluidCpDataTextView.setVisibility(View.GONE);
            fluidPrandtlDataTextView.setVisibility(View.GONE);
            fluidConductivityDataTextView.setVisibility(View.GONE);
            fluidKinematicViscosityDataTextView.setVisibility(View.GONE);
            fluidDynamicViscosityDataTextView.setVisibility(View.GONE);
            fluidDynamicViscositySurfaceTempDataTextView.setVisibility(View.GONE);

            meanTempLabelTextView.setVisibility(View.GONE);
            fluidDensityLabelTextView.setVisibility(View.GONE);
            fluidCpLabelTextView.setVisibility(View.GONE);
            fluidPrandtlLabelTextView.setVisibility(View.GONE);
            fluidConductivityLabelTextView.setVisibility(View.GONE);
            fluidKinematicViscosityLabelTextView.setVisibility(View.GONE);
            fluidDynamicViscosityLabelTextView.setVisibility(View.GONE);
            fluidDynamicViscositySurfaceTempLabelTextView.setVisibility(View.GONE);

        }
        else
        {
            meanTempDataTextView.setVisibility(View.VISIBLE);
            fluidDensityDataTextView.setVisibility(View.VISIBLE);
            fluidCpDataTextView.setVisibility(View.VISIBLE);
            fluidPrandtlDataTextView.setVisibility(View.VISIBLE);
            fluidConductivityDataTextView.setVisibility(View.VISIBLE);
            fluidKinematicViscosityDataTextView.setVisibility(View.VISIBLE);
            fluidDynamicViscosityDataTextView.setVisibility(View.VISIBLE);
            fluidDynamicViscositySurfaceTempDataTextView.setVisibility(View.VISIBLE);

            meanTempLabelTextView.setVisibility(View.VISIBLE);
            fluidDensityLabelTextView.setVisibility(View.VISIBLE);
            fluidCpLabelTextView.setVisibility(View.VISIBLE);
            fluidPrandtlLabelTextView.setVisibility(View.VISIBLE);
            fluidConductivityLabelTextView.setVisibility(View.VISIBLE);
            fluidKinematicViscosityLabelTextView.setVisibility(View.VISIBLE);
            fluidDynamicViscosityLabelTextView.setVisibility(View.VISIBLE);
            fluidDynamicViscositySurfaceTempLabelTextView.setVisibility(View.VISIBLE);

        }


    }

    private void displayLicensesAlertDialog() {
        WebView view = (WebView) LayoutInflater.from(this).inflate(R.layout.dialog_licenses, null);
        view.loadUrl("file:///android_asset/open_source_licenses.html");
        AlertDialog mAlertDialog = new AlertDialog.Builder(this, R.style.Theme_AppCompat_Light_Dialog_Alert)
                .setTitle(getString(R.string.action_licenses))
                .setView(view)
                .setPositiveButton(android.R.string.ok, null)
                .show();
    }


    public void checkLimits(double prefMinValue, double prefMaxValue, String maxLimitString)
    {
        if (prefMaxValue < prefMinValue)
        {
            SharedPreferences sharedPref = PreferenceManager.getDefaultSharedPreferences(this);
            SharedPreferences.Editor editor = sharedPref.edit();
            auxvariable3 = prefMinValue+1;
            editor.putString(maxLimitString,  String.valueOf(auxvariable3));
            editor.commit();
            Toast.makeText(this, R.string.errorMessageLimits, Toast.LENGTH_SHORT).show();
        }
    }

    public void updateLimits()
    {

        SharedPreferences sharedPref = PreferenceManager.getDefaultSharedPreferences(this);

        //updating operating pressure
        inletPressure = Double.parseDouble(sharedPref.getString("inletOperatingPressure", "101.3"))*1000;

        minInletTemperature = Double.parseDouble(sharedPref.getString("minInletTemperature", "1"));
        maxInletTemperature = Double.parseDouble(sharedPref.getString("maxInletTemperature", "100"));
        checkLimits(minInletTemperature, maxInletTemperature, "maxInletTemperature");

        minPipeLength = Double.parseDouble(sharedPref.getString("minPipeLength", "0.01"));
        maxPipeLength = Double.parseDouble(sharedPref.getString("maxPipeLength", "10"));
        checkLimits(minPipeLength, maxPipeLength, "maxPipeLength");

        minSide = Double.parseDouble(sharedPref.getString("minSide", "0.01"));
        maxSide = Double.parseDouble(sharedPref.getString("maxSide", "10"));
        checkLimits(minSide, maxSide, "maxSide");

        minMassFlow = Double.parseDouble(sharedPref.getString("minMassFlow", "0.01"));
        minMassFlow = Double.parseDouble(sharedPref.getString("minMassFlow", "0.01"));
        maxMassFlow = Double.parseDouble(sharedPref.getString("maxMassFlow", "10"));
        checkLimits(minMassFlow, maxMassFlow, "maxMassFlow");

        minVolumetricFlow = Double.parseDouble(sharedPref.getString("minVolumetricFlow", "0.01"));
        maxVolumetricFlow = Double.parseDouble(sharedPref.getString("maxVolumetricFlow", "10"));
        checkLimits(minVolumetricFlow, maxVolumetricFlow, "maxVolumetricFlow");

        minInletVelocity = Double.parseDouble(sharedPref.getString("minInletVelocity", "0.01"));
        maxInletVelocity = Double.parseDouble(sharedPref.getString("maxInletVelocity", "10"));
        checkLimits(minInletVelocity, maxInletVelocity, "maxInletVelocity");

        minSurfaceTemperature = Double.parseDouble(sharedPref.getString("minSurfaceTemperature", "0.01"));
        maxSurfaceTemperature = Double.parseDouble(sharedPref.getString("maxSurfaceTemperature", "10"));
        checkLimits(minSurfaceTemperature, maxSurfaceTemperature, "maxSurfaceTemperature");

        minHeatTransferRate = Double.parseDouble(sharedPref.getString("minHeatTransferRate", "0.01"));
        maxHeatTransferRate = Double.parseDouble(sharedPref.getString("maxHeatTransferRate", "10"));
        checkLimits(minHeatTransferRate, maxHeatTransferRate, "maxHeatTransferRate");


        SeekBar tEntSeekBar = (SeekBar) findViewById(R.id.tempInletSeekBar);
        auxvariable2 = maxInletTemperature - minInletTemperature;
        tEntSeekBar.setMax((int) Math.round(auxvariable2*100));

        SeekBar tubeLSeekBar = (SeekBar) findViewById(R.id.tubeLengthSeekBar);
        auxvariable2 = maxPipeLength - minPipeLength;
        tubeLSeekBar.setMax((int) Math.round(auxvariable2*100));

        SeekBar sideASeekBar = (SeekBar) findViewById(R.id.sideASeekBar);
        SeekBar sideBSeekBar = (SeekBar) findViewById(R.id.sideBSeekBar);
        auxvariable2 = maxSide - minSide;
        sideASeekBar.setMax((int) Math.round(auxvariable2*100));
        sideBSeekBar.setMax((int) Math.round(auxvariable2*100));

        SeekBar flowInputSeekBar = (SeekBar) findViewById(R.id.flowInputSeekBar);
        if(flowInputKnownData == "massFlow")
        {
            auxvariable2 = maxMassFlow - minMassFlow;
        }
        if(flowInputKnownData == "volumetricFlow")
        {
            auxvariable2 = maxVolumetricFlow - minVolumetricFlow;
        }
        if(flowInputKnownData == "meanVelocity")
        {
            auxvariable2 = maxInletVelocity - minInletVelocity;
        }
        flowInputSeekBar.setMax((int) Math.round(auxvariable2*100));

        SeekBar thermalInputSeekBar = (SeekBar) findViewById(R.id.thermalInputSeekBar);
        if(thermalInputKnownData == "heatDissipated")
        {
            auxvariable2 = maxHeatTransferRate - minHeatTransferRate;
        }
        else //surface temperature
        {
            auxvariable2 = maxSurfaceTemperature - minSurfaceTemperature;
        }
        thermalInputSeekBar.setMax((int) Math.round(auxvariable2*100));

        setPropertiesVisibility();

    }


    public void updateCalculation() {

        checkDiameters();
        updateLimits();
        int j=0;
        for(int i=1; i<3; i++)
        {

            //-----Mean temperature calculation----------

            meanTemp = (inletTemp + outletTemp) / 2;

            //------Properties calculation----------------

            fluidDensity = CoolProp.PropsSI("D", "P", inletPressure, "T", meanTemp + 273.15, selectedFluid);//Density en kg/m3
            fluidCp = CoolProp.PropsSI("C", "P", inletPressure, "T", meanTemp + 273.15, selectedFluid);//Cp en J/(kg-K)
            fluidPrandtl = CoolProp.PropsSI("PRANDTL", "P", inletPressure, "T", meanTemp + 273.15, selectedFluid);//Pr
            fluidConductivity = CoolProp.PropsSI("CONDUCTIVITY", "P", inletPressure, "T", meanTemp + 273.15, selectedFluid);//k en W/(m-K)
            fluidDynamicViscosity = CoolProp.PropsSI("V", "P", inletPressure, "T", meanTemp + 273.15, selectedFluid);//Viscosidad dinamica en Pa-s o kg/(m-s)
            fluidDynamicViscositySurfaceTemp = CoolProp.PropsSI("V", "P", inletPressure, "T", surfaceTemp + 273.15, selectedFluid);//Viscosidad dinamica en Pa-s o kg/(m-s)
            fluidKinematicViscosity = fluidDynamicViscosity / fluidDensity;//m2/s

            //--------surfaceArea, sectionArea and Dh Calculation-------------
            if (tubeGeometry == "Circular")
            {
                surfaceArea = Math.PI * Diameter * tubeLength;
                sectionArea = Math.PI * Diameter * Diameter / 4;
                hydraulicDiameter = Diameter;
            }
            if (tubeGeometry == "Rectangular")
            {
                surfaceArea = 2 * tubeLength * (SideA + SideB);
                sectionArea = SideA * SideB;
                hydraulicDiameter = 2 * SideA * SideB / (SideA + SideB);
                if (SideA >= SideB)
                {
                    ratioAB = SideA / SideB;
                } else {
                    ratioAB = SideB / SideA;
                }
            }
            if (tubeGeometry == "Annular")
            {
                surfaceArea = Math.PI * innerDiameter * tubeLength;//Assumes no heat losses to environment
                sectionArea = Math.PI * Math.abs(Diameter * Diameter - innerDiameter*innerDiameter)/ 4;
                hydraulicDiameter =Math.abs(Diameter - innerDiameter);
                ratioDiameters = innerDiameter/Diameter;
            }


            //----------mass flow and Reynolds calculation---------------

            if (flowInputKnownData == "massFlow")
            {
                inletVelocity = flowInputKnownDataValue / (fluidDensity * sectionArea);
            }

            if (flowInputKnownData == "volumetricFlow")
            {
                inletVelocity = flowInputKnownDataValue / sectionArea;
            }

            if (flowInputKnownData == "meanVelocity")
            {
                inletVelocity = flowInputKnownDataValue;
            }

            massFlow = fluidDensity * inletVelocity * sectionArea;
            volumetricFlow = inletVelocity * sectionArea;
            Reynolds = hydraulicDiameter * inletVelocity / fluidKinematicViscosity;

            //------ Flow regime calculation---------

            if (Reynolds < 2300)
            {
                flowRegime = "Laminar";
            }
            else
            {

                if (Reynolds < 10000)
                {
                    flowRegime = "Transition";
                }
                else
                {
                    flowRegime = "Turbulent";
                }

            }

            //------- friction factor calculation ---------

            if (flowRegime == "Laminar")
            {
                if (tubeGeometry == "Rectangular")
                {
                    frictionFactor = (56.05434 * Math.pow(ratioAB, 0.185232)) / Reynolds;
                }
                else //circular or annular
                {
                    frictionFactor = 64 / Reynolds;
                }
            }
            else //turbulent or transition
            {
                if (((tubeMaterial == "smoothTube") || (tubeMaterial == "glass")) || (tubeMaterial == "plastic"))
                {
                    frictionFactor = Math.pow(0.79 * Math.log(Reynolds) - 1.64, -2);
                }
                else
                {
                    frictionFactor = Math.pow(1 / (-1.8 * Math.log10(6.9 / Reynolds + Math.pow(roughness / (3.7 * hydraulicDiameter), 1.11))), 2);
                }
            }


            //----------Entry length calculation---------------

            if (flowRegime == "Laminar")
            {
                hydrodynamicEntryLength = 0.05 * Reynolds * hydraulicDiameter;
                thermalEntryLength = fluidPrandtl * hydrodynamicEntryLength;
            }
            else
            {
                hydrodynamicEntryLength = 10 * hydraulicDiameter;
                thermalEntryLength = hydrodynamicEntryLength;
            }


            //----------Nusselt calculation---------------

            if (flowRegime == "Laminar")
            {
                if (thermalInputKnownData == "heatDissipated")
                {
                    if (tubeGeometry == "Circular")
                    {
                        Nusselt = 4.36;//assumes fully developed flow
                    }

                    if (tubeGeometry == "Rectangular")
                    {
                        Nusselt = -0.0024014931414 * Math.pow(ratioAB, 3) - 0.0071180781406 * Math.pow(ratioAB, 2) + 0.6540471286299 * ratioAB + 2.9344425653701;
                    }

                    if (tubeGeometry == "Annular")
                    {
                        Nusselt = 4.44377320159635 * Math.pow(ratioDiameters, -0.430463997);
                    }

                }
                else//thermalInputKnownData == "surfaceTemperature"
                {
                    if (tubeLength > thermalEntryLength)//Thermally developed flow
                    {
                        if (tubeGeometry == "Circular")
                        {
                            Nusselt = 3.66;
                        }

                        if (tubeGeometry == "Rectangular")
                        {
                            Nusselt = -0.0031930893428 * Math.pow(ratioAB, 3) + 0.0150934478123 * Math.pow(ratioAB, 2) + 0.4742116921787 * ratioAB + 2.468441524652;
                        }

                        if (tubeGeometry == "Annular")
                        {
                            Nusselt = 4.44377320159635 * Math.pow(ratioDiameters, -0.430463997);
                        }

                    }
                    else//thermally entry zone
                    {
                        if (surfaceTemp > (meanTemp + 20))
                        {
                            Nusselt = 1.86 * Math.pow(Reynolds * fluidPrandtl * hydraulicDiameter / tubeLength, 1 / 3) * Math.pow(fluidDynamicViscosity / fluidDynamicViscositySurfaceTemp, 0.14);
                        }
                        else
                        {
                            auxvariable = (hydraulicDiameter / tubeLength) * Reynolds * fluidPrandtl;
                            Nusselt = 3.66 + (0.065 * auxvariable) / (1 + 0.04 * Math.pow(auxvariable, 2 / 3));
                        }
                    }
                }
            }
            else//turbulent or transition
            {
                Nusselt = (frictionFactor / 8) * (Reynolds - 1000) * fluidPrandtl / (1 + 12.7 * Math.pow(frictionFactor / 8, 0.5) * (Math.pow(fluidPrandtl, 2 / 3) - 1));
            }

            //----------Convection coefficient calculation---------------

            convectionCoefficient = Nusselt * fluidConductivity / hydraulicDiameter;

            //----------Outlet temperature, NTU, heat transfer calculation---------------

            NTU = convectionCoefficient * surfaceArea / (massFlow * fluidCp);

            if (thermalInputKnownData == "heatDissipated")
            {
                outletTemp = inletTemp + heatTransferRate / (massFlow * fluidCp);
            }
            else//surfaceTemperature
            {
                outletTemp = surfaceTemp - (surfaceTemp - inletTemp) * Math.exp(-NTU);
                deltaTlm = (inletTemp-outletTemp)/Math.log((surfaceTemp-outletTemp)/(surfaceTemp-inletTemp));
                heatTransferRate =convectionCoefficient*surfaceArea*deltaTlm;
            }
            j=j+1;

        }//End for loop

        //Results update

        outletTempDataTextView.setText(tempFormat.format(outletTemp) + " °C");
        massFlowDataTextView.setText(tempFormat.format(massFlow) + " kg/s");
        volumetricFlowDataTextView.setText(frictionFormat.format(volumetricFlow) + " m³/s");
        inletVelocityDataTextView.setText(frictionFormat.format(inletVelocity) + " m³/s");
        ReynoldsDataTextView.setText(tempFormat.format(Reynolds));
        flowRegimeDataTextView.setText(flowRegime);
        frictionFactorDataTextView.setText(frictionFormat.format(frictionFactor));
        thermalEntryLengthDataTextView.setText(frictionFormat.format(thermalEntryLength)+ " m");
        hydrodynamicEntryLengthDataTextView.setText(frictionFormat.format(hydrodynamicEntryLength)+ " m");
        NusseltDataTextView.setText(tempFormat.format(Nusselt));
        convectionCoefficientDataTextView.setText(tempFormat.format(convectionCoefficient)+ " W/(m²°C)");
        heatTransferRateDataTextView.setText(tempFormat.format(heatTransferRate)+ " W");
        NTUDataTextView.setText(tempFormat.format(NTU));
        meanTempDataTextView.setText(tempFormat.format(meanTemp)+ " °C");
        fluidDensityDataTextView.setText(tempFormat.format(fluidDensity)+ " kg/m³");
        fluidCpDataTextView.setText(tempFormat.format(fluidCp)+ " J/(kg°C)");
        fluidPrandtlDataTextView.setText(tempFormat.format(fluidPrandtl));
        fluidConductivityDataTextView.setText(tempFormat.format(fluidConductivity)+ " W/(m°C)");
        fluidKinematicViscosityDataTextView.setText(viscosityFormat.format(fluidKinematicViscosity)+ " m²/s");
        fluidDynamicViscosityDataTextView.setText(viscosityFormat.format(fluidDynamicViscosity)+ " kg/(m s)");
        fluidDynamicViscositySurfaceTempDataTextView.setText(viscosityFormat.format(fluidDynamicViscositySurfaceTemp)+ " kg/(m s)");
    }


    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);
        Toolbar toolbar = (Toolbar) findViewById(R.id.toolbar);
        setSupportActionBar(toolbar);


        tempFormat.setMinimumFractionDigits(0);
        tempFormat.setMaximumFractionDigits(2);

        lengthFormat.setMinimumFractionDigits(0);
        lengthFormat.setMaximumFractionDigits(2);

        frictionFormat.setMinimumFractionDigits(0);
        frictionFormat.setMaximumFractionDigits(4);

        viscosityFormat.setMinimumFractionDigits(0);
        viscosityFormat.setMaximumFractionDigits(7);

        //Spinner fluido
        Spinner spinnerFluid = (Spinner) findViewById(R.id.spinnerFluido);
        ArrayAdapter<CharSequence> adapter = ArrayAdapter.createFromResource(this, R.array.Fluids, android.R.layout.simple_spinner_item);
        adapter.setDropDownViewResource(android.R.layout.simple_spinner_dropdown_item);
        spinnerFluid.setAdapter(adapter);

        //Spinner tipo de tuberia
        Spinner spinnerGeom = (Spinner) findViewById(R.id.GeometrySpinner);
        ArrayAdapter<CharSequence> adapter1 = ArrayAdapter.createFromResource(this, R.array.Profile, android.R.layout.simple_spinner_item);
        adapter1.setDropDownViewResource(android.R.layout.simple_spinner_dropdown_item);
        spinnerGeom.setAdapter(adapter1);

        //Spinner tube diameter D
        final Spinner spinnerDiam = (Spinner) findViewById(R.id.DiameterSpinner);
        ArrayAdapter<CharSequence> adapter2 = ArrayAdapter.createFromResource(this, R.array.Diameter, android.R.layout.simple_spinner_item);
        adapter2.setDropDownViewResource(android.R.layout.simple_spinner_dropdown_item);
        spinnerDiam.setAdapter(adapter2);

        //Spinner tube diameter d (inner diameter)
        final Spinner spinnerInnerDiam = (Spinner) findViewById(R.id.innerDiameterSpinner);
        ArrayAdapter<CharSequence> adapter21 = ArrayAdapter.createFromResource(this, R.array.innerDiameter, android.R.layout.simple_spinner_item);
        adapter21.setDropDownViewResource(android.R.layout.simple_spinner_dropdown_item);
        spinnerInnerDiam.setAdapter(adapter21);


        //spinnerMaterial
        Spinner spinnerTubeMaterial = (Spinner) findViewById(R.id.spinnerMaterial);
        ArrayAdapter<CharSequence> adapter3 = ArrayAdapter.createFromResource(this, R.array.MaterialList, android.R.layout.simple_spinner_item);
        adapter3.setDropDownViewResource(android.R.layout.simple_spinner_dropdown_item);
        spinnerTubeMaterial.setAdapter(adapter3);

        //spinnerFlowInput
        Spinner spinnerFlowInputInfo = (Spinner) findViewById(R.id.spinnerFlowInput);
        ArrayAdapter<CharSequence> adapter4 = ArrayAdapter.createFromResource(this, R.array.flowInputInfo, android.R.layout.simple_spinner_item);
        adapter4.setDropDownViewResource(android.R.layout.simple_spinner_dropdown_item);
        spinnerFlowInputInfo.setAdapter(adapter4);

        //spinnerThermalInput
        Spinner spinnerThermalInputInfo = (Spinner) findViewById(R.id.spinnerThermalInput);
        ArrayAdapter<CharSequence> adapter5 = ArrayAdapter.createFromResource(this, R.array.thermalInputInfo, android.R.layout.simple_spinner_item);
        adapter5.setDropDownViewResource(android.R.layout.simple_spinner_dropdown_item);
        spinnerThermalInputInfo.setAdapter(adapter5);


        //flowInputDataSeekBar
        flowInputDataTextView = (TextView) findViewById(R.id.flowInputTextView);
        SeekBar flowInputDataSeekBar = (SeekBar) findViewById(R.id.flowInputSeekBar);
        flowInputDataSeekBar.setOnSeekBarChangeListener(flowInputDataSeekBarListener);

        //thermalInputDataSeekBar
        thermalInputDataTextView = (TextView) findViewById(R.id.thermalInputTextView);
        SeekBar thermalInputDataSeekBar = (SeekBar) findViewById(R.id.thermalInputSeekBar);
        thermalInputDataSeekBar.setOnSeekBarChangeListener(thermalInputDataSeekBarListener);

        //tEntSeekBar
        tempinletTextView = (TextView) findViewById(R.id.tempinletTextView);
        SeekBar tEntSeekBar = (SeekBar) findViewById(R.id.tempInletSeekBar);
        tEntSeekBar.setOnSeekBarChangeListener(TempInletSeekBarListener);

        //tubeLSeekBar
        tubeLTextView = (TextView) findViewById(R.id.tubeLengthTextView);
        SeekBar tubeLSeekBar = (SeekBar) findViewById(R.id.tubeLengthSeekBar);
        tubeLSeekBar.setOnSeekBarChangeListener(tubeLengthSeekBarListener);

        //sideASeekBar
        sATextview = (TextView) findViewById(R.id.sideATextview);
        final SeekBar sideA_SeekBar = (SeekBar) findViewById(R.id.sideASeekBar);
        sideA_SeekBar.setOnSeekBarChangeListener(sideASeekBarListener);

        //sideBSeekBar
        sBTextview = (TextView) findViewById(R.id.sideBTextview);
        final SeekBar sideB_SeekBar = (SeekBar) findViewById(R.id.sideBSeekBar);
        sideB_SeekBar.setOnSeekBarChangeListener(sideBSeekBarListener);


        sALabelTextview = (TextView) findViewById(R.id.sideALabelTextview);
        sBLabelTextview = (TextView) findViewById(R.id.sideBLabelTextview);


        //---Results textview declaration-----
        outletTempDataTextView = (TextView) findViewById(R.id.outletTempTextView);
        massFlowDataTextView = (TextView) findViewById(R.id.massFlowTextView);
        volumetricFlowDataTextView = (TextView) findViewById(R.id.volumetricFlowTextView);
        inletVelocityDataTextView = (TextView) findViewById(R.id.inletVelocityTextView);
        ReynoldsDataTextView = (TextView) findViewById(R.id.ReynoldsTextView);
        flowRegimeDataTextView = (TextView) findViewById(R.id.flowRegimeTextView);
        frictionFactorDataTextView = (TextView) findViewById(R.id.frictionFactorTextView);
        thermalEntryLengthDataTextView = (TextView) findViewById(R.id.thermalEntryLengthTextView);
        hydrodynamicEntryLengthDataTextView = (TextView) findViewById(R.id.hydrodynamicEntryLengthTextView);
        NusseltDataTextView = (TextView) findViewById(R.id.NusseltTextView);
        convectionCoefficientDataTextView = (TextView) findViewById(R.id.convectionCoefficientTextView);
        heatTransferRateDataTextView = (TextView) findViewById(R.id.heatTransferRateTextView);
        NTUDataTextView = (TextView) findViewById(R.id.NTUTextView);
        meanTempDataTextView = (TextView) findViewById(R.id.meanTempTextView);
        fluidDensityDataTextView = (TextView) findViewById(R.id.fluidDensityTextView);
        fluidCpDataTextView = (TextView) findViewById(R.id.fluidCpTextView);
        fluidPrandtlDataTextView = (TextView) findViewById(R.id.fluidPrandtlTextView);
        fluidConductivityDataTextView = (TextView) findViewById(R.id.fluidConductivityTextView);
        fluidKinematicViscosityDataTextView = (TextView) findViewById(R.id.fluidKinematicViscosityTextView);
        fluidDynamicViscosityDataTextView = (TextView) findViewById(R.id.fluidDynamicViscosityTextView);
        fluidDynamicViscositySurfaceTempDataTextView = (TextView) findViewById(R.id.fluidDynamicViscositySurfaceTempTextView);


        //---Label results declaration----
        massFlowLabelsTextView = (TextView) findViewById(R.id.massFlowLabelTextView);
        volumetricFlowLabelsTextView = (TextView) findViewById(R.id.volumetricFlowLabelTextView);
        inletVelocityLabelsTextView = (TextView) findViewById(R.id.inletVelocityLabelTextView);

        updateCalculation();

        spinnerTubeMaterial.setOnItemSelectedListener(new AdapterView.OnItemSelectedListener() {

            @Override
            public void onItemSelected(AdapterView<?> adapterView, View view, int position, long id) {

                if (position == 0)
                {
                    tubeMaterial = "smoothTube";
                    roughness = 0;//mm
                }

                if (position == 1)
                {
                    tubeMaterial = "steel";
                    roughness = 0.045;//mm
                }

                if (position == 2)
                {
                    tubeMaterial = "stainlessSteel";
                    roughness = 0.002;//mm
                }

                if (position == 3)
                {
                    tubeMaterial = "glass";
                    roughness = 0;//mm
                }

                if (position == 4)
                {
                    tubeMaterial = "plastic";
                    roughness = 0;//mm
                }

                if (position == 5)
                {
                    tubeMaterial = "wood";
                    roughness = 0.5;//mm
                }

                if (position == 6)
                {
                    tubeMaterial = "rubberSmoothed";
                    roughness = 0.01;//mm
                }

                if (position == 7)
                {
                    tubeMaterial = "copper";
                    roughness = 0.0015;//mm
                }

                if (position == 8)
                {
                    tubeMaterial = "brass";
                    roughness = 0.0015;//mm
                }

                if (position == 9)
                {
                    tubeMaterial = "castIron";
                    roughness = 0.26;//mm
                }

                if (position == 10)
                {
                    tubeMaterial = "galvanizedIron";
                    roughness = 0.15;//mm
                }

                if (position == 11)
                {
                    tubeMaterial = "wrougthIron";
                    roughness = 0.046;//mm
                }

                if (position == 12)
                {
                    tubeMaterial = "concrete";
                    roughness = 1;//mm, range: 0.9 to 9 mm
                }
                roughness =roughness/1000;
                updateCalculation();
            }

            @Override
            public void onNothingSelected(AdapterView<?> parent)
            {
                // vacio
            }
        });//cierra spinnerTubeMaterial


        spinnerFluid.setOnItemSelectedListener(new AdapterView.OnItemSelectedListener() {

            @Override
            public void onItemSelected(AdapterView<?> adapterView, View view, int position, long id) {

                if (position == 0)//Water
                {
                    selectedFluid = "Water";
                }

                if (position == 1)//Aire
                {
                    selectedFluid = "Air";
                }
                updateCalculation();
            }

            @Override
            public void onNothingSelected(AdapterView<?> parent)
            {
                // vacio
            }
        });//cierra spinnerGeom


        spinnerGeom.setOnItemSelectedListener(new AdapterView.OnItemSelectedListener() {

            @Override
            public void onItemSelected(AdapterView<?> adapterView, View view, int position, long id) {
                ImageView img = (ImageView) findViewById(R.id.image);
                diamLabelTextview = (TextView) findViewById(R.id.DiameterLabelTextview);
                innerDiamLabelTextview = (TextView) findViewById(R.id.innerDiameterLabelTextview);
                //Spinner spinnerDiams = (Spinner) findViewById(R.id.DiameterSpinner);
                //Spinner innerSpinnerDiams = (Spinner) findViewById(R.id.innerDiameterSpinner);
                if (position == 0)//Tuberia cilíndrica
                {


                    img.setImageResource(R.drawable.ic_circular);
                    tubeGeometry = "Circular";
                    sALabelTextview.setVisibility(View.GONE);
                    sBLabelTextview.setVisibility(View.GONE);
                    sATextview.setVisibility(View.GONE);
                    sBTextview.setVisibility(View.GONE);
                    sideA_SeekBar.setVisibility(View.GONE);
                    sideB_SeekBar.setVisibility(View.GONE);
                    spinnerDiam.setVisibility(View.VISIBLE);
                    diamLabelTextview.setVisibility(View.VISIBLE);
                    spinnerInnerDiam.setVisibility(View.GONE);
                    innerDiamLabelTextview.setVisibility(View.GONE);

                }

                if (position == 1)//Tuberia rectangular
                {

                    img.setImageResource(R.drawable.ic_rectangular);
                    tubeGeometry = "Rectangular";
                    sALabelTextview.setVisibility(View.VISIBLE);
                    sBLabelTextview.setVisibility(View.VISIBLE);
                    sATextview.setVisibility(View.VISIBLE);
                    sBTextview.setVisibility(View.VISIBLE);
                    sideA_SeekBar.setVisibility(View.VISIBLE);
                    sideB_SeekBar.setVisibility(View.VISIBLE);
                    spinnerDiam.setVisibility(View.GONE);
                    diamLabelTextview.setVisibility(View.GONE);
                    spinnerInnerDiam.setVisibility(View.GONE);
                    innerDiamLabelTextview.setVisibility(View.GONE);
                }

                if (position == 2)//Tuberia anular
                {
                    img.setImageResource(R.drawable.ic_annular);
                    tubeGeometry = "Annular";
                    sALabelTextview.setVisibility(View.GONE);
                    sBLabelTextview.setVisibility(View.GONE);
                    sATextview.setVisibility(View.GONE);
                    sBTextview.setVisibility(View.GONE);
                    sideA_SeekBar.setVisibility(View.GONE);
                    sideB_SeekBar.setVisibility(View.GONE);
                    spinnerDiam.setVisibility(View.VISIBLE);
                    diamLabelTextview.setVisibility(View.VISIBLE);
                    spinnerInnerDiam.setVisibility(View.VISIBLE);
                    innerDiamLabelTextview.setVisibility(View.VISIBLE);

                }
                updateCalculation();
            }

            @Override
            public void onNothingSelected(AdapterView<?> parent)
            {
                // vacio
            }
        });//cierra spinnerGeom

        spinnerDiam.setOnItemSelectedListener(new AdapterView.OnItemSelectedListener() {

            @Override
            public void onItemSelected(AdapterView<?> adapterView, View view, int position, long id) {

                if (position == 0) // 1/8"
                   Diameter = 0.265;

                if (position == 1)// 1/4"
                    Diameter = 0.36;

                if (position == 2)// 3/8"
                    Diameter = 0.489;

                if (position == 3)// 1/2"
                    Diameter = 0.618;

                if (position == 4)// 3/4"
                    Diameter = 0.82;

                if (position == 5)// 1"
                    Diameter = 1.043;

                if (position == 6)// 1-1/4"
                    Diameter = 1.374;

                if (position == 7)// 1-1/2"
                    Diameter = 1.604;

                if (position == 8)// 2"
                    Diameter = 2.059;

                if (position == 9)// 2 1/2"
                    Diameter = 2.459;

                if (position == 10)// 3"
                    Diameter = 3.058;

                if (position == 11)// 3 1/2"
                    Diameter = 3.538;

                if (position == 12)// 4"
                    Diameter = 4.016;

                if (position == 13)// 5"
                    Diameter = 5.037;

                if (position == 14)// 6"
                    Diameter = 6.053;

                if (position == 15)// 8"
                    Diameter = 7.967;

                Diameter = Diameter*25.4/1000;//conversion to m
                updateCalculation();
            }

            @Override
            public void onNothingSelected(AdapterView<?> parent)
            {
                // vacio
            }
        });//cierra spinnerDiam




        spinnerInnerDiam.setOnItemSelectedListener(new AdapterView.OnItemSelectedListener() {

            @Override
            public void onItemSelected(AdapterView<?> adapterView, View view, int position, long id) {

                if (position == 0) // 1/8"
                    innerDiameter = 0.265;

                if (position == 1)// 1/4"
                    innerDiameter = 0.36;

                if (position == 2)// 3/8"
                    innerDiameter = 0.489;

                if (position == 3)// 1/2"
                    innerDiameter = 0.618;

                if (position == 4)// 3/4"
                    innerDiameter = 0.82;

                if (position == 5)// 1"
                    innerDiameter = 1.043;

                if (position == 6)// 1-1/4"
                    innerDiameter = 1.374;

                if (position == 7)// 1-1/2"
                    innerDiameter = 1.604;

                if (position == 8)// 2"
                    innerDiameter = 2.059;

                if (position == 9)// 2 1/2"
                    innerDiameter = 2.459;

                if (position == 10)// 3"
                    innerDiameter = 3.058;

                if (position == 11)// 3 1/2"
                    innerDiameter = 3.538;

                if (position == 12)// 4"
                    innerDiameter = 4.016;

                if (position == 13)// 5"
                    innerDiameter = 5.037;

                if (position == 14)// 6"
                    innerDiameter = 6.053;

                innerDiameter = innerDiameter*25.4/1000;//conversión a metros
                updateCalculation();
            }

            @Override
            public void onNothingSelected(AdapterView<?> parent)
            {
                // vacio
            }
        });//cierra spinnerInnerDiam

        spinnerFlowInputInfo.setOnItemSelectedListener(new AdapterView.OnItemSelectedListener() {

            @Override
            public void onItemSelected(AdapterView<?> adapterView, View view, int position, long id) {

                if(position==0)
                {
                    flowInputKnownData = "massFlow";
                    flowInputDataTextView.setText(tempFormat.format(flowInputKnownDataValue/100) + " kg/s");
                    massFlowLabelsTextView.setVisibility(View.GONE);
                    volumetricFlowLabelsTextView.setVisibility(View.VISIBLE);
                    inletVelocityLabelsTextView.setVisibility(View.VISIBLE);
                    massFlowDataTextView.setVisibility(View.GONE);
                    volumetricFlowDataTextView.setVisibility(View.VISIBLE);
                    inletVelocityDataTextView.setVisibility(View.VISIBLE);
                }

                if(position==1)
                {
                    flowInputKnownData = "volumetricFlow";
                    flowInputDataTextView.setText(tempFormat.format(flowInputKnownDataValue) + " m³/s");
                    massFlowLabelsTextView.setVisibility(View.VISIBLE);
                    volumetricFlowLabelsTextView.setVisibility(View.GONE);
                    inletVelocityLabelsTextView.setVisibility(View.VISIBLE);
                    massFlowDataTextView.setVisibility(View.VISIBLE);
                    volumetricFlowDataTextView.setVisibility(View.GONE);
                    inletVelocityDataTextView.setVisibility(View.VISIBLE);
                }

                if(position==2)
                {
                    flowInputKnownData = "meanVelocity";
                    flowInputDataTextView.setText(tempFormat.format(flowInputKnownDataValue) + " m/s");
                    massFlowLabelsTextView.setVisibility(View.VISIBLE);
                    volumetricFlowLabelsTextView.setVisibility(View.VISIBLE);
                    inletVelocityLabelsTextView.setVisibility(View.GONE);
                    massFlowDataTextView.setVisibility(View.VISIBLE);
                    volumetricFlowDataTextView.setVisibility(View.VISIBLE);
                    inletVelocityDataTextView.setVisibility(View.GONE);
                }
                updateCalculation();
            }

            @Override
            public void onNothingSelected(AdapterView<?> parent)
            {
                // vacio
            }
        });//cierra spinnerFlowInputInfo


        spinnerThermalInputInfo.setOnItemSelectedListener(new AdapterView.OnItemSelectedListener() {

            @Override
            public void onItemSelected(AdapterView<?> adapterView, View view, int position, long id) {

                if(position==0)
                {
                    thermalInputKnownData = "surfaceTemperature";
                    thermalInputDataTextView.setText(tempFormat.format(thermalInputKnownDataValue+20) + " °C");
                    surfaceTemp=thermalInputKnownDataValue+20;

                }

                if(position==1)
                {
                    thermalInputKnownData = "heatDissipated";
                    thermalInputDataTextView.setText(tempFormat.format(thermalInputKnownDataValue) + " W");
                }

                updateCalculation();
            }

            @Override
            public void onNothingSelected(AdapterView<?> parent)
            {
                // vacio
            }
        });//cierra spinnerThermalInputInfo

    }//Cierra onCreate

    @Override
    public boolean onCreateOptionsMenu(Menu menu) {
        // Inflate the menu; this adds items to the action bar if it is present.
        getMenuInflater().inflate(R.menu.menu_main, menu);
        return true;
    }

    @Override
    public void onResume(){
        super.onResume();
        updateCalculation();
    }

    @Override
    public boolean onOptionsItemSelected(MenuItem item) {
        // Handle action bar item clicks here. The action bar will
        // automatically handle clicks on the Home/Up button, so long
        // as you specify a parent activity in AndroidManifest.xml.

        switch (item.getItemId()) {
            case R.id.helpLicense:
                displayLicensesAlertDialog();
                return true;
            default:
                Intent preferencesIntent = new Intent(this, SettingsActivity.class);
                startActivity(preferencesIntent);
                return super.onOptionsItemSelected(item);
        }


    }


    //seekbar de dato conocido de entrada del fluido (flujo másico, caudal o velocidad de entrada)
    private final SeekBar.OnSeekBarChangeListener flowInputDataSeekBarListener =
            new SeekBar.OnSeekBarChangeListener() {

                @Override
                public void onProgressChanged(SeekBar seekBar, int progress, boolean fromUser) {

                    flowInputKnownDataValue = progress;
                    if(flowInputKnownData == "massFlow")
                    {
                        flowInputKnownDataValue = (float)progress/100 + minMassFlow;
                        flowInputDataTextView.setText(tempFormat.format(flowInputKnownDataValue) + " kg/s");
                    }
                    if(flowInputKnownData == "volumetricFlow")
                    {
                        flowInputKnownDataValue = (float)progress/100 + minVolumetricFlow;
                        flowInputDataTextView.setText(tempFormat.format(flowInputKnownDataValue) + " m³/s");

                    }
                    if(flowInputKnownData == "meanVelocity")
                    {
                        flowInputKnownDataValue = (float)progress/100 + minInletVelocity;
                        flowInputDataTextView.setText(tempFormat.format(flowInputKnownDataValue) + " m/s");
                    }
                    updateCalculation();
                }
                @Override
                public void onStartTrackingTouch(SeekBar seekBar) { }
                @Override
                public void onStopTrackingTouch(SeekBar seekBar) { }
            };

    //seekbar de dato conocido térmico de entrada del fluido (temperatura superficial o calor disipado)
    private final SeekBar.OnSeekBarChangeListener thermalInputDataSeekBarListener =
            new SeekBar.OnSeekBarChangeListener() {

                @Override
                public void onProgressChanged(SeekBar seekBar, int progress, boolean fromUser) {

                    if(thermalInputKnownData == "surfaceTemperature")
                    {
                        thermalInputKnownDataValue = (float)progress/100 + minSurfaceTemperature;
                        thermalInputDataTextView.setText(tempFormat.format(thermalInputKnownDataValue) + " °C");
                        surfaceTemp = thermalInputKnownDataValue;

                    }
                    if(thermalInputKnownData == "heatDissipated")
                    {
                        thermalInputKnownDataValue = (float)progress/100 + minHeatTransferRate;
                        thermalInputDataTextView.setText(tempFormat.format(thermalInputKnownDataValue) + " W");
                        heatTransferRate = thermalInputKnownDataValue;
                    }
                    updateCalculation();
                }
                @Override
                public void onStartTrackingTouch(SeekBar seekBar) { }
                @Override
                public void onStopTrackingTouch(SeekBar seekBar) { }
            };


    //seekbar del lado A del ducto rectangular
    private final SeekBar.OnSeekBarChangeListener sideASeekBarListener =
            new SeekBar.OnSeekBarChangeListener() {

                @Override
                public void onProgressChanged(SeekBar seekBar, int progress, boolean fromUser) {
                    SideA = (float)progress/100+0.01;
                    sATextview.setText(lengthFormat.format(SideA) + "m");
                    updateCalculation();

                }
                @Override
                public void onStartTrackingTouch(SeekBar seekBar) { }
                @Override
                public void onStopTrackingTouch(SeekBar seekBar) { }
            };


    //seekbar del lado B del ducto rectangular
    private final SeekBar.OnSeekBarChangeListener sideBSeekBarListener =
            new SeekBar.OnSeekBarChangeListener() {

                @Override
                public void onProgressChanged(SeekBar seekBar, int progress, boolean fromUser) {
                    SideB = (float)progress/100+0.01;
                    sBTextview.setText(lengthFormat.format(SideB) + "m");
                    updateCalculation();

                }
                @Override
                public void onStartTrackingTouch(SeekBar seekBar) { }
                @Override
                public void onStopTrackingTouch(SeekBar seekBar) { }
            };

    //seekbar de temperatura de entrada a la tuberia
    private final SeekBar.OnSeekBarChangeListener TempInletSeekBarListener =
            new SeekBar.OnSeekBarChangeListener() {

                @Override
                public void onProgressChanged(SeekBar seekBar, int progress, boolean fromUser) {
                    inletTemp = (float)progress/100+minInletTemperature;
                    tempinletTextView.setText(tempFormat.format(inletTemp) + "°C");
                    updateCalculation();

                }
                @Override
                public void onStartTrackingTouch(SeekBar seekBar) { }
                @Override
                public void onStopTrackingTouch(SeekBar seekBar) { }
            };

    //seekbar de longitud de tuberia
    private final SeekBar.OnSeekBarChangeListener tubeLengthSeekBarListener =
            new SeekBar.OnSeekBarChangeListener() {

                @Override
                public void onProgressChanged(SeekBar seekBar, int progress, boolean fromUser) {
                    tubeLength = (float)progress/100 + minPipeLength;
                    tubeLTextView.setText(lengthFormat.format(tubeLength) + " m");
                    updateCalculation();
                }
                @Override
                public void onStartTrackingTouch(SeekBar seekBar) { }
                @Override
                public void onStopTrackingTouch(SeekBar seekBar) { }
            };


}//cierra MainActivity
