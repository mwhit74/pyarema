# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 12:51:05 2020

@author: mwhitten
"""

import math

def eq13141(fax, fa_ax, fb1, fa_b1, fb2, fa_b2, k1, l1, r1, 
            k2, l2, r2, e):
    """Axial Compression and Bending
    
    AREMA 2018 Section 1.3.14.1
    
    Member subject to both axial compression and bending stresses shall be
    proportioned to satisfy the following requirements:
        
        
    if f_axial/fa_axial <= 0.15:
        IR = fa/Fa + fb1/Fb1 + fb2/Fb2
    else:
        moment_mag1 = 1-f_axial/(0.514*math.pow(math.pi,2)*e)*math.pow(klr1,2)
        moment_mag2 = 1-f_axial/(0.514*math.pow(math.pi,2)*e)*math.pow(klr2,2)
        
        IR = fa/Fa + fb1/(Fb1*moment_mag1) + fb2/(Fb2*moment_mag2)
        
        
    Args:
        fax (float): calculated axial stress [psi]
        
        fb1 (float): calculated bending stess about axis 1-1 [psi]
        
        fb2 (float): calculated bending stress about axis 2-2 [psi]
        
        fa_ax (flaot): allowable axial stress if axial force alone 
                       existed [psi]
                          
        fa_b1 (float): allowable compression bending stress about axis 1-1,
                       that would be permitted if bending alone existed [psi]
                       
        fa_b2 (float): allowable compression bending stress about axis 2-2,
                       that would be permitted if bending alone existed [psi]
                       
        k1 (float): effective length factor about 1-1 axis
        
        l1 (float): unbraced length about 1-1 axis [inches]
        
        r1 (float): radius of gyration about 1-1 axis [inches]
        
        k2 (float): effective length factor about 2-2 axis
        
        l2 (float): unbraced length about 2-2 axis [inches]
        
        r2 (float): radius of gyration about 2-2 axis [inches]
        
        e (float): modulus of elasticity [psi]
        
    Returns:
        ir (tuple(float, str)): interaction ratio of the axial stress and 
                                bi-axail bending stress per the prescribed 
                                equations
    
    Notes:
        1. Units are in lbs and inches.
        2. To ignore moment magnification factors enter k1 = k2 = 0.0.
           Secondary effects, where they are determined to be significant,
           should be considered elsewhere in the analysis or design if the 
           moment magnification is ignored here. 
    """
    
    ref_text = "AREMA 2018 Section 1.3.14.1 \n\n"
    
    user_input = (f'fax = {fax:.2f}, fa_ax = {fa_ax:.2f},' +
                   f'fb1 = {fb1:.2f}, fa_b1 = {fa_b1:.2f},' +
                   f'fb2 = {fb2:.2f}, fa_b2 = {fa_b2:.2f},' +
                   f'k1 = {k1:.2f}, l1 = {l1:.2f}, r1 = {r1:.2f},' +
                   f'k2 = {k2:.2f}, l2 = {l2:.2f}, r2 = {r2:.2f},' +
                   f'E = {e:.2f} \n\n')
    
    fax_ir = fax/fa_ax
    fb1_ir = fb1/fa_b1
    fb2_ir = fb2/fa_b2
    
    text1 = (f'fa/Fa = {fax:.2f}/{fa_ax:.2f} \n' +
             f'fa/Fa = {fax_ir:.2f} \n\n')
    text2 = (f'fb1/Fb1 = {fb1:.2f}/{fa_b1:.2f} \n' +
             f'fb1/Fb1 = {fb1_ir:.3f} \n\n')
    text3 = (f'fb2/Fb2 = {fb2:.2f}/{fa_b2:.2f} \n' +
             f'fb2/Fb2 = {fb2_ir:.3f} \n\n')
    
    if fax_ir <= 0.15:
        ir = fax_ir + fb1_ir + fb2_ir
        
        text4 = (f'IR = fax/Fax + fb1/Fb1 + fb2/Fb2 \n'
                 f'IR = {fax_ir:.2f} + {fb1_ir:.2f} + {fb2_ir:.2f} \n'
                 f'IR = {ir:.3f}')
        
        text = text1 + text2 + text3 + text4
        
    else: 
        klr1 = k1*l1/r1
        klr2 = k2*l2/r2
    
        moment_mag1 = 1-fax/(0.514*math.pow(math.pi,2)*e)*math.pow(klr1,2)
        moment_mag2 = 1-fax/(0.514*math.pow(math.pi,2)*e)*math.pow(klr2,2)
        
        ir = fax_ir + fb1_ir/moment_mag1 + fb2/moment_mag2
        
        text5 = f'k1*l1/r1 = {klr1:.2f} \n'
        text6 = f'k2*l2/r2 = {klr2:.2f} \n\n'
    
        text7 = (f'Moment Mag. 1 = 1-fax/' +
                 f'(0.514*math.pow(math.pi,2)*e)*math.pow(k1*l1/r1,2) \n' +
                 f'Moment Mag. 1 = 1-{fax:.2f}/' + 
                 f'(0.514*math.pow(math.pi,2)*{e:.1f})'
                 f'*math.pow({klr1:.2f},2) \n' +
                 f'Moment Mag. 1 = {moment_mag1:.3f} \n\n')
        text8 = (f'Moment Mag. 2 = 1-fax/' +
                 f'(0.514*math.pow(math.pi,2)*e)*math.pow(k2*l2/r2,2) \n' +
                 f'Moment Mag. 2 = 1-{fax:.2f}/' + 
                 f'(0.514*math.pow(math.pi,2)*{e:.1f})'
                 f'*math.pow({klr2:.2f},2) \n' +
                 f'Moment Mag. 2 = {moment_mag2:.3f} \n\n')
        
        text9 = (f'IR = fax/Fax + fb1/(Fb1*Moment Mag. 1)'
                 f' + fb2/(Fb2*Moment Mag. 2) \n'
                 f'IR = {fax_ir:.2f} + {fb1_ir/moment_mag1:.2f}'
                 f' + {fb2_ir/moment_mag1:.2f} \n'
                 f'IR = {ir:.3f} \n')
        
        text = text1 + text2 + text3 + text5 + text6 + text7 + text8 + text9
        
    text = ref_text + user_input + text
        
    return ir, text
        

def eq141d1(stress_condition, fy, fu):
    """Axial tension allowable stress subject to direct tensile load
    
    AREMA 2018 Table 15-1-11 Row 1 Equations
    
    
    if stress_condition == 1:
        fa_t = 0.55*fy
    elif stress_condition == 2:
        fa_tn = 0.47*fu
    elif stress_condition == 3:
        fa_tnp = 0.45*fy
    
    
    Args:
        stress_condition (int): 
            1 - Axial tension, structural steel, gross section
            2 - Axial tension, structural steel, effective net section
                (see articles 1.5.8 and 1.6.5)
            3 - Axial tension, structural steel, effective net area
                at cross-section of pin hole of pin-connected members
                
        fy (float): yield stress [psi]
        
        fu (float): ultimate stress [psi]
        
    Returns:
        fa_t (tuple(float, str)): allowable tension stress on gross 
                                  section [psi]
        
        -OR-
        
        fa_tn (tuple(float, str)): allowable tension stress on net 
                                   section [psi]
        
        -OR-
        
        fa_tnp (tuple(float, str)): allowable tension stress on net section at
                                    cross-section of pin hole of pin-connected
                                    members [psi]
                                    
    Notes:
        1. Units in lbs and inches.
    
    """
    
    ref_text = "AREMA 2018 Section 1.4.1 Table 15-1-11 Row 1 \n\n"
    
    user_input = (f'Stress Condition = {stress_condition:d},' +
                  f'Fy = {fy:.2f}, Fu = {fu:.2f} \n')
    
    if stress_condition == 1:
        fa = 0.55*fy
        text = (f'1 - Axial tension, structural steel, gross section \n' +
                f'fa_tg = 0.55*fy \n' +
                f'fa_tg = 0.55*{fy:.1f} \n' +
                f'fa_tg = {fa:.1f}')
    elif stress_condition == 2:
        fa = 0.47*fu
        text = (f'2 - Axial tension, structural steel, effective net section' +
                f'(see articles 1.5.8 and 1.6.5) \n'
                f'fa_tn = 0.47*fu \n' +
                f'fa_tn = 0.47*{fu:.1f} \n' +
                f'fa_tn = {fa:.1f}')
    elif stress_condition == 3:
        fa = 0.45*fy
        text = (f'3 - Axial tension, structural steel, effective net area' +
                f'at cross-section of pin hole of pin-connected members \n'
                f'fa_tnp = 0.45*fy \n' +
                f'fa_tnp = 0.45*{fy:.1f} \n' +
                f'fa_tnp = {fa:.1f}')
        
    text = ref_text + user_input + text
        
    return fa, text

def eq141d2():
    pass

def eq141d3(fy):
    """Tension in extreme fibers subject to bending, net section
    
    AREMA 2018 Table 15-1-11 Row 3 Equations
    
    Tension in extreme fibers of rolled shapes, girders and built-up section, 
    subject to bending, net section
    
    
    fa_bt = 0.55*fy
    
    
    Args:                
        fy (float): yield stress [psi]

        
    Returns:
        fa_bt (tuple(float, str)): allowable bending tensile stress [psi]
    
    """
    
    ref_text = "AREMA 2018 Section 1.4.1 Table 15-1-11 Row 3 \n\n"
    
    user_input = f'Fy = {fy:.2f} \n'
    
    fa_bt = 0.55*fy
    
    text = (f'Tension in extreme fibers of rolled shapes, girder and \n' +
            f'built-up section, subject to bending, net section \n' +
            f'fa_bt = 0.55*fy \n' +
            f'fa_bt = 0.55*{fy:.1f} \n' +
            f'fa_bt = {fa_bt:.1f}')
    
    text = ref_text + user_input + text
    
    return fa_bt, text

def eq141d4():
    pass

def eq141d5(k, l, r, fy, e):
    """Axial compression allowable stress, gross section
    
    AREMA 2018 Table 15-1-11 Row 5 Equations
    (Except for stiffeners of beams and girders, and splice material)
    
    
    klr = k*l/r
    elastic_limit = 0.629*math.sqrt(fy/e)
    plastic_limit = 5.034*math.sqrt(fy/e)
    
    if klr <= elastic_limit:
        fa_ac = 0.55*fy
    elif klr > elastic_limit and klr <= plastic_limit:
        fa_ac = 0.60*fy-math.pow(17500.0*fy/e,3/2)*klr
    elif klr > plastic_limit:
        fa_ac = 0.514*math.pow(math.pi, 2)*e/math.pow(klr,2)
    
    
    Args:
        k (float): effective length factor in the axis considered
                   k = 7/8 for members with pin-end connections
                   k = 3/4 for members with bolted or welded end connections
        
        l (float): unbraced length in the axis considered [inches]
        
        r (float): radius of gyration in the axis considered [inches]
        
        fy (float): yield stress [psi]
        
        e (float): modulus of elasticity [psi]
        
    Returns:
        fa_ac (tuple(float, str)): allowable axial compression on gross section
    
    Notes:
        1. Units are in lbs and inches.
    """
    
    ref_text = "AREMA 2018 Section 1.4.1 Table 15-1-11 Row 5 \n\n"
    
    user_input = (f'k = {k:.2f}, l = {l:.2f}, r = {r:.2f}, Fy = {fy:.2f},' +
                  f'E = {e:.2f} \n')
    
    klr = k*l/r
    elastic_limit = 0.629/math.sqrt(fy/e)
    plastic_limit = 5.034/math.sqrt(fy/e)
    
    if klr <= elastic_limit:
        fa_ac = 0.55*fy
        text = (f'kl/r <= 0.629/math.sqrt(fy/e) \n' +
                f'{klr:.2f} <= {elastic_limit:.2f} \n' + 
                f'fa_ac = 0.55*fy \n' + 
                f'fa_ac = 0.55*{fy:.1f} \n' +
                f'fa_ac = {fa_ac:.1f}')
    elif klr > elastic_limit and klr <= plastic_limit:
        fa_ac = 0.60*fy-math.pow(17500.0*fy/e,3/2)*klr
        text = (f'kl/r > 0.629/math.sqrt(fy/e)'
                f'and kl/r <= 5.034/math.sqrt(fy/e) \n' + 
                f'{klr:.2f} > {elastic_limit:.2f}'
                f'and {klr:.2f} <= {plastic_limit:.2f} \n' +             
                f'fa_ac = 0.60*fy-math.pow(17500.0*fy/e,3/2)*kl/r \n' + 
                f'fa_ac = 0.60*{fy:.1f}-'
                f'math.pow(17500.0*{fy:.1f}/{e:.1f},3/2)*{klr:.2f} \n' +
                f'fa_ac = {fa_ac:.1f}')
    elif klr > plastic_limit:
        fa_ac = 0.514*math.pow(math.pi, 2)*e/math.pow(klr,2)
        text = (f'kl/r > 5.034/math.sqrt(fy/e) \n' +
                f'{klr:.2f} > {plastic_limit:.2f} \n' +
                f'fa_ac = 0.514*math.pow(math.pi, 2)*e/math.pow(kl/r,2) \n' + 
                f'fa_ac = 0.514*math.pow(math.pi, 2)*{e:.2f}' +
                f'/math.pow({klr:.2f},2) \n' +
                f'fa_ac = {fa_ac:.1f}')
        
    text = ref_text + user_input + text
    
    return fa_ac, text
        
def eq141d6(fy):
    """Compression in extreme fibers of I-type section bending about weak axis
    
    AREMA 2018 Section 1.4.1 Table 15-1-11 Row 6
    
    Compression in extreme fibers of I-type members subjected to loading
    perpendicular to the web.
    
    
    fa_b = 0.55*fy
    
    
    Args:
        fy (float): yield stress [psi]
        
    Return:
        fa_bt (tuple(float, str)): allowable bending tensile stress [psi]
    
    """
    
    ref_text = "AREMA 2018 Section 1.4.1 Table 15-1-11 Row 6 \n\n"
    
    user_input = (f'Fy = {fy:.2f} \n\n')


    fa_b = 0.55*fy
    
    
    text = (f'fa_b = 0.55*fy \n' +
            f'fa_b = 0.55*{fy:.1f} \n' +
            f'fa_b = {fa_b:.1f}')
    
    text = ref_text + user_input + text
    
    return fa_b, text
    

def eq141d7(l, ryc, afc, d, fy, e):
    """Compression in extreme fibers of flexural members allowable stress
    
    AREMA 2018 Table 15-1-11 Row 7 Equations
    
    Compression in extreme fibers of flexural members symmertical about the
    principal axis in the plane of the web (other than box-type flexural 
    members) that are rolled beams or welded built-up members with solid
    rectangular flanges, the larger of the values computer by the following
    
    
    fa_bc1 = (0.55*fy - 0.55*math.pow(fy,2)/
           (6.3*math.pow(math.pi,2)*e)*math.pow(l/ryc,2)
           
    fa_bc2 = 0.131*math.pi*e/(l*d*math.sqrt(1+mu)/afc)
    
    fa_bc3 = 0.55*fy
    
    fa_bc = min(max(fa1, fa2), fa3)
    
    
    Args:
        l (float): distance bwteen points of lateral support for the
                   compression flange, unbraced length [inches]
        
        ryc (float): minimum radius of gyration of the compression flange and
                    that portion of the web area on the compression side of the
                    axis of bending, about an axis in the plane of the 
                    web [inches]
                    
        afc (float): area of the smaller flange excluding any portion of the 
                    web [inches^2]
        
        d (float): overall dpeth of the member [inches]
        
        fy (float): yield stress [psi]
        
        e (float): modulus of elasticity [psi]
        
    Returns:
        fa_bc (float): allowable compression stress in extreme fibers of 
                    flexural members [psi]
        
    Notes:
        1. Units are in lbs and inches.
        2. Poisson's ratio, mu, is taken as 0.3.
    
    """
    ref_text = "AREMA 2018 Section 1.4.1 Table 15-1-11 Row 7 \n\n"
    
    user_input = (f'l = {l:.2f}, ryc = {ryc:.2f}, afc = {afc:.2f}, ' +
                  f'd = {d:.2f}, Fy = {fy:.2f}, E = {e:.2f} \n\n')
    
    mu = 0.3
    
    fa_bc1 = (0.55*fy - 0.55*math.pow(fy,2)/
           (6.3*math.pow(math.pi,2)*e)*math.pow(l/ryc,2))
           
    fa_bc2 = 0.131*math.pi*e/(l*d*math.sqrt(1+mu)/afc)
    
    fa_bc3 = 0.55*fy
    
    fa_bc = min(max(fa_bc1, fa_bc2), fa_bc3)
    
    text1 = (f'fa_bc1 = 0.55*fy - 0.55*math.pow(fy,2)/' +
             f'(6.3*math.pow(math.pi,2)*e)*math.pow(l/ryc,2) \n' +
             f'fa_bc1 = 0.55*{fy:.1f} - 0.55*math.pow({fy:.1f},2)/' +
             f'(6.3*math.pow(math.pi,2)*{e:.1f})' +
             f'*math.pow({l:.2f}/{ryc:.2f},2) \n' +
             f"""fa_bc1 = {fa_bc1:.1f} \n\n""")
    
    text2 = (f'fa_bc2 = 0.131*math.pi*e/(l*d*math.sqrt(1+mu)/afc) \n' +
             f'fa_bc2 = 0.131*math.pi*{e:.1f}/' +
             f'({l:.2f}*{d:.2f}*math.sqrt(1+{mu:.2f})/{afc:.2f}) \n' +
             f'fa_bc2 = {fa_bc2:.1f} \n\n')
    
    text3 = (f'fa_bc3 = 0.55*fy \n' +
             f'fa_bc3 = 0.55*{fy:.2f} \n' +
             f'fa_bc3 = {fa_bc3:.1f} \n\n')
    
    text4 = (f'fa_bc = min(max(fa_bc1, fa_bc2), fa_bc3) \n' +
             f'fa_bc = min(max({fa_bc1:.2f}, {fa_bc2:.2f}), {fa_bc3:.2f}) \n' +
             f'fa_bc = min({max(fa_bc1, fa_bc2):.2f}, {fa_bc3:.2f}) \n' +
             f'fa_bc = {min(max(fa_bc1, fa_bc2), fa_bc3):.1f} \n')
    
    text = ref_text + user_input + text1  + text2 + text3 + text4
    
    return fa_bc, text

def eq141d10(l, sx, a, sum_st, iy, fy, e):
    """Compression in extreme fibers of box type flexural members
    
    AREMA 2018 Section 1.4.1 Table 15-1-11 Row 10
    
    Compression in the extreme fibers of box type welded or bolted flexural
    members symmetrical about the principal axis midway between the webs
    
    
    (l/r)e = math.sqrt(1.105*math.pi/sxx*sqrt(sum(s/t))/
             a*math.sqrt(i_yy/(1+mu)))
    
    fa_bc = 0.55*fy - 0.55*math.pow(fy,2)/(6.3*math.pow(math.pi,2)*e)*
            math.pow((l/r)e,2)
            
            
    Args:
        l (float): distance between points of lateral support for the
                   compression flange, unbraced length [inches]
                  
        sx (float): section modulus of the box type member about its
                    major axis [inches^3]
                    
        a (float): total area enclosed within the center lines of the box
                   type member webs and flanges [inches^2]
                   
        sum_st (float): sum of the ratio width-to-thickness of each flange and 
                        ratio of the depth to thickness of each web (neglect
                        any portion of the flange which projects beyond the 
                        box section)
                        
        iy (float): second moment of area of the  box type member about its
                    minor axis, [inches^4]
            
    Returns:
        fa_bc (float): allowable compression stress in extreme fibers of
                       box type flexure members
    
    Notes:
        1. Units in lbs and inches.
        2. Poisson's ratio, mu, is taken as 0.3.
    """
    
    ref_text = "AREMA 2018 Section 1.4.1 Table 15-1-11 Row 10 \n\n"
    
    user_input = (f'l = {l:.2f}, Sx = {sx:.2f}, a = {a:.2f}, ' + 
                 f'sum_st = {sum_st:.2f}, Iy = {iy:.2f}, Fy = {fy:.1f}, ' +
                 f'E = {e:.1f} \n\n')
    
    mu = 0.3

    lre = math.sqrt((1.105*math.pi/sx*math.sqrt(sum_st))/
                    (a*math.sqrt(iy/(1+mu))))
    
    fa_bc = (0.55*fy-0.55*math.pow(fy,2)/
             (6.3*math.pow(math.pi,2)*e)*math.pow(lre,2))
    
    text1 = (f'(l/r)e = math.sqrt((1.105*math.pi/sx*math.sqrt(sum_st))/' +
            f'(a*math.sqrt(iy/(1+mu)))) \n' +
            f'(l/r)e = math.sqrt((1.105*math.pi/{sx:.2f}*math.sqrt({sum_st:.2f}))/' +
            f'({a:.2f}*math.sqrt({iy:.2f}/(1+{mu:.2f})))) \n' +
            f'(l/r)e = {lre:.2f} \n')
            
    text2 = (f'fa_bc = (0.55*fy-0.55*math.pow(fy,2)/' +
            f'(6.3*math.pow(math.pi,2)*e)*math.pow(lre,2)) \n' +
            f'fa_bc = (0.55*{fy:.1f}-0.55*math.pow({fy:.1f},2)/' +
            f'(6.3*math.pow(math.pi,2)*{e:.1f})*math.pow({lre:.2f},2)) \n' +
            f'fa_bc = {fa_bc:.1f}')
            
    text = ref_text + user_input + text1 + text2
            
    return fa_bc, text

def eq141d13(fy):
    """Shear in webs of rolled beams and plate girders, gross section
    
    AREMA 2018 Table 15-1-11 Row 13 Equations 
    
    
    fa = 0.35*fy
    
    
    Args:
        fy (float): yield stress
        
    Returns:
        fa (tuple(float, str)): Allowable shear stress in webs
    
    """
    ref_text = "AREMA 2018 Section 1.4.1 Table 15-1-11 Row 13 \n\n"
    
    user_input = f'fy = {fy:.1f} psi \n\n' 
    
    fa_v = 0.35*fy
    
    text = (f'fa_v = 0.35*fy \n' +
            f'fa_v = 0.35*{fy:.2f} \n' +
            f'fa_v = {fa_v:.2f} \n')
    
    text = ref_text + user_input + text
    
    return fa_v, text

def eq141d17(l, d, fu):
    """Bearing on F3125 Grade A325 and Grade A490 bolts
    
    AREMA 2018 Section 1.4.1 Table 15-1-11 Row 17
    
    Allowable stress of bearing material on F3125 Gr A325 and Gr A490 bolts.
    
    
    fa_bolt_brg = min(l*fu/(2*d), 1.2*fu)
    
    
    Args:
        l (float): distance measured in the line of the force form the center 
                   line of a bolt to the nearest edge of an adjacent bolt or to
                   the end of the connected part toward which the force is
                   directed [inches]
        
        d (float): diameter of bolt [inches]
        
        fu (float): lowest specified minimum tensile strength of the
                    connected part [psi]
                    
    Returns:
        fa_bolt_brg (tuple(float, str)): allowable bearing stress 
        
    """
    ref_text = 'AREMA 2018 Section 1.4.1 Table 15-1-11 Row 17  \n\n'
    
    user_input = f'l = {l:.2f}, d = {d:.2f}, fu = {fu:.1f} \n\n'
    
    fa_brg_bolt = min(l*fu/(2*d), 1.2*fu)
    
    text = (f'fa_brg_bolt = min(l*fu/(2*d), 1.2*fu) \n' +
            f'fa_brg_bolt = min({l:.2f}*{fu:.1f}/(2*{d:.2f}),' 
            f'1.2*{fu:.1f}) \n'+
            f'fa_brg_bolt = min({l*fu/(2*d):.1f}, {1.2*fu:.1f}) \n' +
            f'fa_brg_bolt = {fa_brg_bolt:.1f}')
    
    text = ref_text + user_input + text

    return fa_brg_bolt, text
    

def tb15111a(hole_type, surface_type, bolt_grade):
    """Allowable stress for slip-critical connections
    
    Allowable stress for slip-critical connections
    (Slip Load per Unit Area of Bolt)
    
    Args:
        hole_type (int):
            1 - standard and short slotted perpendicular to the direction
                of load
            2 - Oversize and short slotted parallel to the direction of load
            3 - Long slotted any direction
        surface_type (str):
            A - Class A (slip coefficient 0.30) Uncoated clean mill scale and
                blast-cleaned surfaces with Class A coatings
            B - Class B (slip coefficient 0.50) Uncoated blast-cleaned surfaces
                with Class B coatings and unsealed thermal-sprayed surfaces
            C - Class C (slip coefficient 0.30) Hot-dip galvanized
            D - Class D (slip coefficient 0.45) Blast-cleaned surfaces with 
                Class D coatings
        bolt_grade (int):
            1 - F3125 Grade A325
            2 - F3125 Grade A490
            
    Returns:
        fa_bolt (tuple(float, str)): allowable stress for slip 
                                   critical connection [psi]
    
    """
    
    fa_code = str(bolt_grade) + surface_type + str(hole_type)
    
    fa_bolt = {'1A1':12900,'1B1':21500,'1C1':12900,'1D1':19300,
               '2A1':11000,'2B1':18300,'2C1':11000,'2D1':16500,
               '3A1':9000,'3B1':15000,'3C1':9000,'3D1':13500,
               '1A2':16200,'1B2':26900,'1C2':16200,'1D2':24200,
               '2A2':13800,'2B2':23000,'2C2':13800,'2D2':20700,
               '3A2':11300,'3B2':18900,'3C2':11300,'3D2':17000,}
    
    ht_text = {1:('1 - standard and short slotted perpendicular to the '+
                 'direction of load'),
               2:('2 - Oversize and short slotted parallel to the direction ' +
                 'of load'),
               3:'3 - Long slotted any direction'}
    
def eq162a(element_type, case, t, fy, e):
    """Outstanding element in compression
    
    The width of outstanding elements of members in compression shall not 
    exceed the following, where t, inches, is the thickness of the element:
        
    """
    pass
    