from model import *
from utilities.is800_2007 import IS800_2007
from utilities.other_standards import IS1363_part_1_2002, IS1363_part_3_2002, IS1367_Part3_2002
from utilities.common_calculation import *
import math
import logging
flag = 1
logger = None
beam_d = 0
beam_B = 0

def module_setup():
    global logger
    logger = logging.getLogger("osdag.Tension_calc")


module_setup()

#######################################################################

# cl 6.2 Design Strength Due to Yielding of Gross Section
def tension_member_design_due_to_yielding_of_gross_section(A_g, F_y):
    "design strength of members under axial tension,T_dg,as governed by yielding of gross section"
    "A_g = gross area of cross-section"
    "gamma_m0 = partial safety factor for failure in tension by yielding"
    "F_y = yield stress of the material"
    gamma_m0 = IS800_2007.cl_5_4_1_Table_5["gamma_m0"]['yielding']
    T_dg = A_g * F_y / gamma_m0

    return T_dg



# cl 6.2 Design Strength Due to Rupture of critical section
def preliminary_tension_member_design_due_to_rupture_of_critical_section(A_n, F_u, no_of_bolts):
    "preliminary design strength,T_pdn,as governed by rupture at net section"
    "A_n = net area of the total cross-section"
    "A_nc = net area of the connected leg"
    "A_go = gross area of the outstanding leg"
    "alpha_b,alpha_w = 0.6 - two bolts, 0.7 - three bolts or 0.8 - four or more bolts/welded"
    "gamma_m1 = partial safety factor for failure in tension by ultimate stress"
    "F_u = Ultimate Strength of material"
    "w = outstanding leg width"
    "b_s = shear lag width"
    "t = thickness of the leg"
    "Lc = length of the end connection"

    # if connection_type == "bolted":
    #     if no_of_bolts == "2":
    #         alpha = 0.6
    #     elif no_of_bolts == "3":
    #         alpha = 0.7
    #     else:
    #         aplha = 0.8
    # else:
    #     alpha = 0.8

    if no_of_bolts <=2:
        alpha = 0.6
    elif no_of_bolts == 3:
        alpha = 0.7
    else:
        alpha = 0.8

    gamma_m1 = IS800_2007.cl_5_4_1_Table_5["gamma_m1"]['ultimate_stress']
    T_pdn = alpha * A_n * F_u / gamma_m1

    return T_pdn

# cl 6.2 Design Strength Due to Rupture of critical section
def tension_member_design_due_to_rupture_of_critical_section(A_n, F_u):
    "preliminary design strength,T_pdn,as governed by rupture at net section"
    "A_n = net area of the total cross-section"
    "A_nc = net area of the connected leg"
    "A_go = gross area of the outstanding leg"
    "alpha_b,alpha_w = 0.6 - two bolts, 0.7 - three bolts or 0.8 - four or more bolts/welded"
    "gamma_m1 = partial safety factor for failure in tension by ultimate stress"
    "F_u = Ultimate Strength of material"
    "w = outstanding leg width"
    "b_s = shear lag width"
    "t = thickness of the leg"
    "Lc = length of the end connection"

    # if connection_type == "bolted":
    #     if no_of_bolts == "2":
    #         alpha = 0.6
    #     elif no_of_bolts == "3":
    #         alpha = 0.7
    #     else:
    #         aplha = 0.8
    # else:
    #     alpha = 0.8

    gamma_m1 = IS800_2007.cl_5_4_1_Table_5["gamma_m1"]['ultimate_stress']
    T_pdn = 0.9 * A_n * F_u / gamma_m1

    return T_pdn

def tension_angle_member_design_due_to_rupture_of_critical_section(A_nc, A_go, F_u, F_y, L_c, w, b_s, t):
    "design strength,T_dn,as governed by rupture at net section"
    "A_n = net area of the total cross-section"
    "A_nc = net area of the connected leg"
    "A_go = gross area of the outstanding leg"
    "alpha_b,alpha_w = 0.6 - two bolts, 0.7 - three bolts or 0.8 - four or more bolts/welded"
    "gamma_m1 = partial safety factor for failure in tension by ultimate stress"
    "F_u = Ultimate Strength of material"
    "w = outstanding leg width"
    "b_s = shear lag width"
    "t = thickness of the leg"
    "L_c = length of the end connection"
    "gamma_m0 = partial safety factor for failure in tension by yielding"
    "F_y = yield stress of the material"

    gamma_m0 = IS800_2007.cl_5_4_1_Table_5["gamma_m0"]['yielding']
    gamma_m1 = IS800_2007.cl_5_4_1_Table_5["gamma_m1"]['ultimate_stress']

    beta = float(1.4 - (0.076 * float(w) / float(t) * float(F_y) / float(F_u) * float(b_s) / float(L_c)))
    print(beta)

    if beta <= (F_u * gamma_m0 / F_y * gamma_m1) and beta >= 0.7:
        beta = beta
    else:
        beta = 0.7

    T_dn = (0.9 * A_nc * F_u / gamma_m1) + (beta * A_go * F_y / gamma_m0)

    return T_dn

def tension_member_design_check_for_slenderness(K,L,r,load_type="none"):
    "A tension member in which a reversal of direct stress occurs due to loads otherthan wind or seismic forces less than 180"
    "KL= effective length of member"
    "r = radius of gyration of member"

    # if load_type == "Reversal Load":
    #     if K * L / r < 180:
    #         design_check = True
    #     else:
    #         design_check = False
    # else:
    #     if K * L / r < 400:
    #         design_check = True
    #     else:
    #         design_check = False
    slender = K * L*1000 /r

    return slender

# Start of Main Program

def tension_design(uiObj):
    global logger
    global design_status
    design_status = True

    # if uiObj['Member']['Location'] == "Web":
    #     conn_type = "Web"
    # elif uiObj['Member']['Location'] == "Flange":
    #     conn_type = "Flange"
    # else:
    #     conn_type = "Leg"
    conn = uiObj['Member']['Location']
    Member_type = uiObj['Member']['SectionType']
    # Member_type = "Angles"
    Member_size = uiObj['Member']['SectionSize']
    # Member_size = "40 40 x 4"
    Member_fu = float(uiObj['Member']['fu (MPa)'])
    Member_fy = float(uiObj['Member']['fy (MPa)'])
    Member_length = float(uiObj["Member"]["Member_length"])
    Tension_load = float(uiObj["Load"]["AxialForce (kN)"])

    bolt_dia = int(uiObj['Bolt']['Diameter (mm)'])
    bolt_row = int(uiObj["Bolt"]["RowsofBolts"])
    bolt_column = int(uiObj["Bolt"]["ColumnsofBolts"])
    bolt_row_pitch = float(uiObj["Bolt"]["Rowpitch"])
    bolt_column_pitch = float(uiObj["Bolt"]["Columnpitch"])
    bolt_enddistance = float(uiObj["Bolt"]["Enddistance"])
    bolt_edgedistance = float(uiObj["Bolt"]["Edgedistance"])

    dia_hole = bolt_dia + int(uiObj["bolt"]["bolt_hole_clrnce"])
    end1_cond1 = uiObj["Support_Condition"]["end1_cond1"]
    end1_cond2 = uiObj["Support_Condition"]["end1_cond2"]
    end2_cond1 = uiObj["Support_Condition"]["end2_cond1"]
    end2_cond2 = uiObj["Support_Condition"]["end2_cond2"]

    # mu_f = float(uiObj["bolt"]["slip_factor"])
    # gamma_mw = float(uiObj["weld"]["safety_factor"])
    #
    # gamma_mw = float(uiObj["weld"]["safety_factor"])
    # if gamma_mw == 1.50:
    #     weld_fabrication = 'field'
    # else:
    #     weld_fabrication = 'shop'
    #
    # dp_bolt_hole_type = uiObj["bolt"]["bolt_hole_type"]
    # if dp_bolt_hole_type == "Over-sized":
    #     bolt_hole_type = 'over_size'
    # else:  # "Standard"
    #     bolt_hole_type = 'standard'
    #
    # dia_hole = bolt_dia + int(uiObj["bolt"]["bolt_hole_clrnce"])
    #
    # if uiObj["detailing"]["typeof_edge"] == "a - Sheared or hand flame cut":
    #     edge_type = 'hand_flame_cut'
    # else:   # "b - Rolled, machine-flame cut, sawn and planed"
    #     edge_type = 'machine_flame_cut'
    #
    # corrosive_influences = False
    # if uiObj['detailing']['is_env_corrosive'] == "Yes":
    #     corrosive_influences = True
    #
    # [bolt_shank_area, bolt_net_area] = IS1367_Part3_2002.bolt_area(bolt_dia)

    old_beam_section = get_oldbeamcombolist()
    old_column_section = get_oldcolumncombolist()

    if Member_size in old_beam_section or Member_size in old_column_section:
        logger.warning(": You are using a section (in red colour) that is not available in the latest version of IS 808")

    if Member_fu < 410 or Member_fy < 230:
        logger.warning(" : You are using a section of grade that is not available in the latest version of IS 2062")

    dictmemberdata = get_memberdata(Member_size,Member_type)
    print dictmemberdata
    if Member_type != "Angles":
        member_tw = float(dictmemberdata["tw"])
        member_tf = float(dictmemberdata["T"])
        member_d = float(dictmemberdata["D"])
        member_B = float(dictmemberdata["B"])
        member_R1= float(dictmemberdata["R1"])
        Member_Ag = float (dictmemberdata["Area"]) * 100
        radius_gyration = min((float(dictmemberdata["rz"])),(float(dictmemberdata["ry"])))*10
    else:
        member_leg = dictmemberdata["AXB"]
        leg = member_leg.split("x")
        leg1 = leg[0]
        leg2 = leg[1]
        t = float(dictmemberdata["t"])
        Member_Ag = float(dictmemberdata["Area"]) * 100
        radius_gyration = min((float(dictmemberdata["ru(max)"])), (float(dictmemberdata["rv(min)"]))) * 10

    # Calculation for Design Strength Due to Yielding of Gross Section
    tension_yielding = tension_member_design_due_to_yielding_of_gross_section(Member_Ag,Member_fy)/1000
    no_of_bolts = bolt_row

    k =  IS800_2007.effective_length_coefficeint(end1_cond1, end1_cond2, end2_cond1, end2_cond2)
    tension_slenderness = tension_member_design_check_for_slenderness(k,Member_length,radius_gyration)
    radius_gyration_min = k * (Member_length * 1000)/400

    if conn == "Web":
        Member_An = Member_Ag - (bolt_column * dia_hole * member_tw)
        if bolt_column >=2:
            A_vg = ((bolt_row_pitch*(bolt_row-1) + bolt_enddistance)*member_tw)*2
            A_vn = ((bolt_row_pitch*(bolt_row-1) + bolt_enddistance - ((bolt_row -0.5)*dia_hole)) * member_tw)*2
            A_tg = bolt_column_pitch* (bolt_column - 1)  * member_tw
            A_tn = (bolt_column_pitch*(bolt_column-1) - ((1)*dia_hole)) * member_tw
        else:
            A_vg = (bolt_row_pitch + bolt_enddistance) * member_tw
            A_vn = ((bolt_row_pitch * (bolt_row - 1) + bolt_enddistance - ((bolt_row - 0.5) * dia_hole)) * member_tw)
            A_tg = (bolt_column_pitch + bolt_edgedistance)* member_tw
            A_tn = (bolt_column_pitch + bolt_edgedistance - 0.5* dia_hole) * member_tw
    elif conn == "Back to Back Web":
        Member_Ag = float (dictmemberdata["Area"]) * 100 *2
        Member_An = Member_Ag - (bolt_column * dia_hole * member_tw)
        if bolt_column >=2:
            A_vg = ((bolt_row_pitch*(bolt_row-1) + bolt_enddistance)*member_tw)*2*2
            A_vn = ((bolt_row_pitch*(bolt_row-1) + bolt_enddistance - ((bolt_row -0.5)*dia_hole)) * member_tw)*2*2
            A_tg = bolt_column_pitch* (bolt_column - 1)  * member_tw*2
            A_tn = (bolt_column_pitch*(bolt_column-1) - ((1)*dia_hole)) * member_tw*2
        else:
            A_vg = (bolt_row_pitch + bolt_enddistance) * member_tw*2
            A_vn = ((bolt_row_pitch * (bolt_row - 1) + bolt_enddistance - ((bolt_row - 0.5) * dia_hole)) * member_tw)*2
            A_tg = (bolt_column_pitch + bolt_edgedistance)* member_tw*2
            A_tn = (bolt_column_pitch + bolt_edgedistance - 0.5* dia_hole) * member_tw*2
    elif conn=="Flange":
        if Member_type == "Beams":
            Member_An = Member_Ag - (bolt_column * dia_hole * member_tf)
            A_vg = ((bolt_row_pitch * (bolt_row - 1) + bolt_enddistance) * member_tf) * 2 * 2
            A_vn = ((bolt_row_pitch * (bolt_row - 1) + bolt_enddistance - ((bolt_row - 0.5) * dia_hole)) * member_tf) * 2 * 2
            A_tg = (bolt_column_pitch*(bolt_column/2 -1) + bolt_edgedistance) * member_tf * 2 * 2
            A_tn = ((bolt_column_pitch*(bolt_column/2 -1) + bolt_edgedistance) - (bolt_column/2 -0.5) * dia_hole) * member_tf * 2 * 2
        elif Member_type == "Channels":
            Member_An = Member_Ag - (bolt_column * dia_hole * member_tf)
            A_vg = ((bolt_row_pitch * (bolt_row - 1) + bolt_enddistance) * member_tf) * 2
            A_vn = ((bolt_row_pitch * (bolt_row - 1) + bolt_enddistance - ((bolt_row - 0.5) * dia_hole)) * member_tf) * 2
            A_tg = (bolt_column_pitch * (bolt_column / 2 - 1) + bolt_edgedistance) * member_tf * 2
            A_tn = ((bolt_column_pitch * (bolt_column / 2 - 1) + bolt_edgedistance) - (bolt_column / 2 - 0.5) * dia_hole) * member_tf * 2
    elif conn == "Back to Back Leg":
        Member_Ag = float(dictmemberdata["Area"]) * 100*2
        Member_An = Member_Ag - (bolt_column * dia_hole * t)
        A_vg = ((bolt_row_pitch * (bolt_row - 1) + bolt_enddistance) * t)*2
        A_vn = ((bolt_row_pitch * (bolt_row - 1) + bolt_enddistance - ((bolt_row - 0.5) * dia_hole)) * t)*2
        A_tg = ((bolt_column_pitch * (bolt_column - 1)) + bolt_edgedistance)* t*2
        A_tn = ((((bolt_column_pitch * (bolt_column - 1)) + bolt_edgedistance)) - (((bolt_column -0.5) * dia_hole))) * t*2
    else:
        Member_An = Member_Ag - (bolt_column * dia_hole * t)
        A_vg = ((bolt_row_pitch * (bolt_row - 1) + bolt_enddistance) * t)
        A_vn = ((bolt_row_pitch * (bolt_row - 1) + bolt_enddistance - ((bolt_row - 0.5) * dia_hole)) * t)
        A_tg = ((bolt_column_pitch * (bolt_column - 1)) + bolt_edgedistance)* t
        A_tn = ((((bolt_column_pitch * (bolt_column - 1)) + bolt_edgedistance)) - (((bolt_column -0.5) * dia_hole))) * t

    if Member_type == "Angles":
        w = max(float(leg1), float(leg2))
        shear_lag = ((min(float(leg1), float(leg2)))-bolt_edgedistance) + w - t
        L_c = (bolt_row_pitch * (bolt_row - 1))
    else:
        pass

    # cl 6.2 Design Strength Due to Block Shear
    tension_blockshear = IS800_2007.cl_6_4_1_block_shear_strength(A_vg, A_vn, A_tg, A_tn, Member_fu, Member_fy)/1000
    # Calculation for Design Strength Due to Yielding of Gross Section
    if (Member_type == "Angles")and bolt_row > 1 :
        tension_rupture = tension_angle_member_design_due_to_rupture_of_critical_section(Member_An,Member_Ag, Member_fu , Member_fy, L_c, w, shear_lag, t)/1000
    else:
        tension_rupture = tension_member_design_due_to_rupture_of_critical_section(Member_An, Member_fu)/1000

    tension_design = min(tension_blockshear,tension_rupture,tension_yielding)


 # End of Calculation, SAMPLE Output dictionary
    outputobj = dict()

    # FOR OUTPUT DOCK
    outputobj['Tension_Force'] = {}

    outputobj['Tension_Force']['Yielding'] = float(round(tension_yielding,3))
    outputobj['Tension_Force']['Rupture'] = float(round(tension_rupture,3))
    outputobj['Tension_Force']['Block_Shear'] = float(round(tension_blockshear,3))
    outputobj['Tension_Force']['Efficiency'] = float(round((Tension_load/tension_design),3))
    outputobj['Tension_Force']['Slenderness'] = float(round((tension_slenderness),3))

    # for i,j in uiObj.items():
    #     if j == " ":
    #         logger.error(": Please enter all the inputs")
    #     else:
    #         pass

    if outputobj['Tension_Force']['Efficiency'] < 1 and outputobj['Tension_Force']['Slenderness'] < 400:
        design_status = True
    elif outputobj['Tension_Force']['Efficiency'] > 1 and outputobj['Tension_Force']['Slenderness'] < 400:
        design_status = False
        logger.error(": Chosen Member Section Size is not sufficient")
        logger.info(": Increase the size of Member ")
    elif outputobj['Tension_Force']['Efficiency'] < 1 and outputobj['Tension_Force']['Slenderness']> 400:
        design_status = False
        logger.error(": Chosen Member Section Size is not sufficient")
        logger.warning(": Minimum Radius of Gyration of Member shall be {} mm ".format(radius_gyration_min))
        logger.info(": Increase the size of Member ")

    if design_status is True:
        logger.info(":Member is safe for the applied tension load \n")
        logger.info(":In case of reversal load, Slenderness value shall be less than 180 \n")
        logger.debug(" :=========End Of design===========")
    else:
        logger.error(":Member fails for the applied tension load \n ")
        logger.debug(" :=========End Of design===========")

    return outputobj

# tg = tension_member_design_due_to_yielding_of_gross_section(586, 210)
# print(tg)
#
# tr = preliminary_tension_member_design_due_to_rupture_of_critical_section(586, 410, "welded", "none")
# print(tr)
#
# trc = tension_member_design_due_to_rupture_of_critical_section(300, 300, 410, 210, 50, 50, 50, 6)
# print(trc)
#
# tdb = IS800_2007.cl_6_4_1_block_shear_strength(300, 300, 0, 0, 410, 210)
# print(tdb)
#
# tc = tension_member_design_check_for_slenderness(1,2000,21.1,"Only_tension")
# print(tc)