import math

# --- C++ helper function equivalents ---
def myasinh(X):
    return math.log(X + math.sqrt(X*X + 1.0))

def myacosh(X):
    # Ensure X >= 1 to avoid domain errors, allowing for small precision issues
    if X < 1.0:
        if X > 1.0 - 1e-9: # If very close to 1, treat as 1
            X = 1.0
        else:
            raise ValueError(f"Domain error for myacosh: X = {X} < 1")
    return math.log(X + math.sqrt(X*X - 1.0))

# --- Constants and Input Parameters ---
FP_COORDS = {"x": 5.2, "y": 0.0, "z": -70.0}
AP_COORDS = {"x": 853.87, "y": 0.0, "z": -320.0}
L_TOTAL = 902.2  # m
# RHO_LINE = 697.69  # kg/m (Weight per unit length w = RHO_LINE * g is implicitly handled if Ans_x gives H/wh)
XACC = 1e-4      # Calculation accuracy

MAXIT_P0_SOLVER = 100 # Max iterations for p0 Newton loop in funcd_p0_solver
MAXIT_RTSAFE = 200    # Max iterations for rtsafe

# --- Geometric Parameters from Inputs ---
# h_span is 'h' in the C++ code context
h_span = abs(FP_COORDS["z"] - AP_COORDS["z"])
# L_APFP is the horizontal distance between Fairlead and Anchor
L_APFP = math.sqrt((FP_COORDS["x"] - AP_COORDS["x"])**2 +
                   (FP_COORDS["y"] - AP_COORDS["y"])**2)

# Parameters 'd_param' and 'l_param' for funcd, as per C++ ModuleCatenary::AssRes
# L_total is the full length of the line.
d_param = (L_APFP - (L_TOTAL - h_span)) / h_span
l_param = L_TOTAL / h_span

# --- Core C++ Logic Replication: funcd and rtsafe ---

# This global variable will hold the p0 value associated with the current x_val in funcd,
# similar to how a class member might behave in C++.
# It's primarily for ensuring the internal logic for p0 in funcd is correctly structured.
_current_p0_for_funcd = 0.0

def _solve_for_p0_in_funcd(x_candidate, l_val, xacc_p0):
    """Solves for p0 using Newton's method, internal to funcd logic."""
    p0 = 0.0  # Initial guess from C++ code's behavior
    for _ in range(MAXIT_P0_SOLVER):
        cos_p0 = math.cos(p0)
        if abs(cos_p0) < 1e-12: # Avoid division by zero if cos_p0 is extremely small
            # This implies p0 is near +/- pi/2, which can make terms undefined.
            # Depending on context, could return last p0 or signal an issue.
            # print(f"Warning: cos(p0) is near zero ({cos_p0}) in _solve_for_p0_in_funcd for x={x_candidate}.")
            break 

        tan_p0 = math.tan(p0)
        func1 = 1.0 / x_candidate + 1.0 / cos_p0
        
        sqrt_arg_f1 = func1**2 - 1.0
        # Robust sqrt: handle small negative due to precision, ensure non-negative for real sqrt
        if sqrt_arg_f1 < 0 and sqrt_arg_f1 > -1e-9: # If slightly negative due to precision
            sqrt_arg_f1 = 0.0
        elif sqrt_arg_f1 < 0:
            # print(f"Warning: Negative sqrt_arg_f1 ({sqrt_arg_f1}) in p0 solver for x={x_candidate}, p0={p0}. Returning current p0.")
            return p0 # Cannot proceed with real sqrt

        sqrt_val_f1 = math.sqrt(sqrt_arg_f1)
        f1_for_p0 = x_candidate * (sqrt_val_f1 - tan_p0) - l_val

        if abs(f1_for_p0) < xacc_p0:
            return p0

        # Derivative df1 for p0 Newton step
        if abs(sqrt_val_f1) < 1e-9: # Avoid division by zero for derivative term
            # If sqrt_val_f1 is zero, implies func1 is +/-1. Derivative might be infinite.
            # A robust Newton solver might switch to bisection or use a damped step.
            # For now, if this happens and f1_for_p0 isn't small, we can't improve p0.
            # print(f"Warning: sqrt_val_f1 is near zero in p0 solver derivative for x={x_candidate}, p0={p0}.")
            break

        term_A_df1 = func1 * tan_p0 / (cos_p0 * sqrt_val_f1)
        term_B_df1 = 1.0 / (cos_p0**2) # since tan(p0)^2 + 1 = 1/cos(p0)^2
        df1_for_p0 = x_candidate * (term_A_df1 - term_B_df1)

        if abs(df1_for_p0) < 1e-12: # Denominator for Newton step too small
            # print(f"Warning: df1_for_p0 is near zero ({df1_for_p0}) in p0 solver for x={x_candidate}, p0={p0}.")
            break
        
        p0 = p0 - f1_for_p0 / df1_for_p0
        
        # Keep p0 in a typical range like (-pi/2, pi/2) to avoid issues with tan/cos
        if not (-math.pi/2 * 0.99 < p0 < math.pi/2 * 0.99) :
            # print(f"Warning: p0 ({p0}) moved to edge of range in p0 solver. Clamping or stopping.")
            break # Or clamp p0, though C++ doesn't show clamping.
            
    # print(f"p0 solver for x={x_candidate} finished. p0={p0}, f1_for_p0={f1_for_p0 if 'f1_for_p0' in locals() else 'N/A'}")
    return p0


def _funcd_equations(x_val, d_val, l_val, xacc_for_p0_solver):
    """Calculates f(x) and df(x) based on C++ funcd logic."""
    global _current_p0_for_funcd # To simulate class member update
    f_res, df_res = 0.0, 0.0
    p0_local = 0.0

    if abs(x_val) < 1e-9: # x_val is effectively zero
        f_res = -d_val
        df_res = 0.0 # As in C++
        p0_local = 0.0
    elif x_val > 0:
        if l_val <= 0: # This case is unlikely for typical mooring lines (L_TOTAL > 0, h_span > 0)
            X1 = 1.0 / x_val + 1.0
            term_acosh_X1 = myacosh(X1) # myacosh handles X1 slightly < 1
            sqrt_1_2x = math.sqrt(1.0 + 2.0 * x_val)
            f_res = x_val * term_acosh_X1 - sqrt_1_2x + 1.0 - d_val
            
            sqrt_X1sq_m1 = math.sqrt(max(0,X1**2 - 1.0))
            if abs(x_val * sqrt_X1sq_m1) < 1e-9: # Avoid division by zero in df
                df_res = term_acosh_X1 - (1.0 / sqrt_1_2x) # Simplified derivative
            else:
                df_res = term_acosh_X1 - (1.0 / sqrt_1_2x) - (1.0 / (x_val * sqrt_X1sq_m1))
            p0_local = 0.0
        else: # l_val > 0
            threshold_x = (l_val**2 - 1.0) / 2.0
            if x_val > threshold_x:
                p0_local = _solve_for_p0_in_funcd(x_val, l_val, xacc_for_p0_solver)
                tan_p0 = math.tan(p0_local)
                X2 = l_val / x_val + tan_p0
                asinh_X2 = myasinh(X2)
                asinh_tan_p0 = myasinh(tan_p0)
                f_res = x_val * (asinh_X2 - asinh_tan_p0) - l_val + 1.0 - d_val
                
                sqrt_X2sq_p1 = math.sqrt(X2**2 + 1.0)
                if abs(x_val * sqrt_X2sq_p1) < 1e-9:
                    df_res = asinh_X2 - asinh_tan_p0
                else:
                    df_res = (asinh_X2 - asinh_tan_p0) - l_val / (x_val * sqrt_X2sq_p1)
            else: # x_val <= threshold_x
                X5 = 1.0 / x_val + 1.0
                term_acosh_X5 = myacosh(X5)
                sqrt_1_2x_alt = math.sqrt(1.0 + 2.0 * x_val)
                f_res = x_val * term_acosh_X5 - sqrt_1_2x_alt + 1.0 - d_val

                sqrt_X5sq_m1 = math.sqrt(max(0,X5**2 - 1.0))
                if abs(x_val * sqrt_X5sq_m1) < 1e-9:
                    df_res = term_acosh_X5 - (1.0 / sqrt_1_2x_alt)
                else:
                    df_res = term_acosh_X5 - (1.0 / sqrt_1_2x_alt) - (1.0 / (x_val * sqrt_X5sq_m1))
                p0_local = 0.0
    else: # x_val < 0
        raise ValueError("x_val < 0 encountered in _funcd_equations")
    
    _current_p0_for_funcd = p0_local # Update the 'member-like' p0
    return f_res, df_res, p0_local


def _rtsafe_solver(x1_rt, x2_rt, xacc_rt, d_rt, l_rt, xacc_p0_rt):
    """Solves for root using rtsafe algorithm (Newton-Raphson with bisection backup)."""
    p0_for_root = 0.0 # Stores p0 associated with the final root
    
    f_low, _, p0_low = _funcd_equations(x1_rt, d_rt, l_rt, xacc_p0_rt)
    f_high, _, p0_high = _funcd_equations(x2_rt, d_rt, l_rt, xacc_p0_rt)

    if (f_low > 0 and f_high > 0) or (f_low < 0 and f_high < 0):
        raise ValueError(f"Root not bracketed in rtsafe: f({x1_rt})={f_low}, f({x2_rt})={f_high}")

    if f_low == 0: p0_for_root = p0_low; return x1_rt, p0_for_root
    if f_high == 0: p0_for_root = p0_high; return x2_rt, p0_for_root

    if f_low < 0:
        xl_rt, xh_rt = x1_rt, x2_rt
    else:
        xl_rt, xh_rt = x2_rt, x1_rt
    
    rts_curr = 0.5 * (x1_rt + x2_rt)
    dx_old = abs(x2_rt - x1_rt)
    dx_curr = dx_old
    
    f_curr, df_curr, p0_curr = _funcd_equations(rts_curr, d_rt, l_rt, xacc_p0_rt)
    p0_for_root = p0_curr

    for _ in range(MAXIT_RTSAFE):
        if (((rts_curr - xh_rt) * df_curr - f_curr) * ((rts_curr - xl_rt) * df_curr - f_curr) > 0.0) or \
           (abs(2.0 * f_curr) > abs(dx_old * df_curr)):
            # Bisection step
            dx_old = dx_curr
            dx_curr = 0.5 * (xh_rt - xl_rt)
            rts_curr = xl_rt + dx_curr
            if xl_rt == rts_curr: # Root found (interval small enough)
                # Ensure p0_for_root is for this final rts_curr
                _, _, p0_for_root = _funcd_equations(rts_curr, d_rt, l_rt, xacc_p0_rt)
                return rts_curr, p0_for_root
        else:
            # Newton-Raphson step
            dx_old = dx_curr
            if abs(df_curr) < 1e-12: # Derivative too small, avoid division by zero
                # Fallback to bisection-like step if Newton is unstable
                dx_curr = 0.5 * (xh_rt - xl_rt)
                rts_curr = xl_rt + dx_curr
                if xl_rt == rts_curr: 
                     _, _, p0_for_root = _funcd_equations(rts_curr, d_rt, l_rt, xacc_p0_rt)
                     return rts_curr, p0_for_root
            else:
                dx_curr = f_curr / df_curr

            temp_rts = rts_curr
            rts_curr -= dx_curr
            if temp_rts == rts_curr: # No change, root found
                _, _, p0_for_root = _funcd_equations(rts_curr, d_rt, l_rt, xacc_p0_rt)
                return rts_curr, p0_for_root
        
        if abs(dx_curr) < xacc_rt:
            _, _, p0_for_root = _funcd_equations(rts_curr, d_rt, l_rt, xacc_p0_rt)
            return rts_curr, p0_for_root

        f_curr, df_curr, p0_curr = _funcd_equations(rts_curr, d_rt, l_rt, xacc_p0_rt)
        p0_for_root = p0_curr # Keep track of p0 for current rts_curr

        if f_curr < 0:
            xl_rt = rts_curr
        else:
            xh_rt = rts_curr
            
    raise RuntimeError("Max iterations exceeded in rtsafe.")

# --- Main Calculation Logic ---
print("--- mooring line initial shape calculation ---")
print(f"Fairlead: {FP_COORDS}")
print(f"Anchor: {AP_COORDS}")
print(f"Line Total Length (L_total): {L_TOTAL} m")
print(f"Vertical Span (h_span): {h_span:.3f} m")
print(f"Horizontal Span (L_APFP): {L_APFP:.3f} m")
print(f"Parameter d: {d_param:.5f}")
print(f"Parameter l: {l_param:.5f}")

ans_x_final = 0.0
# p0_final_at_root = 0.0 # Not strictly needed for coordinates but part of C++ state

# Check C++ edge cases before rtsafe
if d_param <= 0:
    ans_x_final = 0.0
    # _current_p0_for_funcd = 0.0 (already default or set by last funcd call)
    print("Condition d_param <= 0 met. Ans_x set to 0.")
elif l_param >= 1.0 and d_param >= (math.sqrt(l_param**2 - 1.0) - (l_param - 1.0)):
    # This specific error check implies l_param >= 1 for sqrt(l_param^2-1)
    raise ValueError("Line length L_total is too short for the separation (C++ error condition).")
else:
    # Initial bracket for rtsafe. Test f(small_x) and f(large_x).
    # Based on manual tests, f(0.001) was neg, f(100) was pos for the example.
    rt_x1 = XACC # Small positive to avoid x=0 issues if any
    rt_x2 = max(100.0, l_param * l_param * 2.0) # Heuristic for upper bound

    f1_check, _, _ = _funcd_equations(rt_x1, d_param, l_param, XACC)
    f2_check, _, _ = _funcd_equations(rt_x2, d_param, l_param, XACC)

    if (f1_check * f2_check) > 0:
        print(f"Warning: Initial rtsafe bracket [{rt_x1}, {rt_x2}] gives f signs: ({f1_check}, {f2_check}). Trying wider.")
        rt_x1 = XACC
        rt_x2 = 1.0e6 # Wide range from C++ example
        f1_check_wide, _, _ = _funcd_equations(rt_x1, d_param, l_param, XACC)
        f2_check_wide, _, _ = _funcd_equations(rt_x2, d_param, l_param, XACC)
        if (f1_check_wide * f2_check_wide) > 0:
             raise ValueError("Cannot bracket root for rtsafe even with wide range. Problem with inputs or funcd.")
        else:
             print(f"Using wide bracket for rtsafe: [{rt_x1}, {rt_x2}]")
    
    print(f"Solving for Ans_x with rtsafe: bracket [{rt_x1:.3e}, {rt_x2:.3e}], xacc={XACC}")
    ans_x_final, _ = _rtsafe_solver(rt_x1, rt_x2, XACC, d_param, l_param, XACC)
    # p0_final_at_root = p0_val_from_rtsafe # Store if needed later

print(f"Calculated Ans_x: {ans_x_final:.6f}")

# Catenary parameter 'a' = H/w = Ans_x * h_span
a_catenary = ans_x_final * h_span
print(f"Catenary parameter 'a': {a_catenary:.3f} m")

# S_suspended_by_ans_x = h_span * math.sqrt(1 + 2 * ans_x_final)
# print(f"Suspended length implied by Ans_x and h_span: {S_suspended_by_ans_x:.3f} m")
# print(f"(Compare with L_TOTAL = {L_TOTAL} m)")
# For node distribution, we use L_TOTAL as the arc length to segment.

# --- Determine global catenary geometry (lowest point x_gm, z_offset) ---
# Line from FP (xf, zf) to AP (xa, za)
xf, zf = FP_COORDS["x"], FP_COORDS["z"]
xa, za = AP_COORDS["x"], AP_COORDS["z"]

# Signed distances: AP relative to FP
Xs_rel = xa - xf  # Horizontal component
Zs_rel = za - zf  # Vertical component
S_line = L_TOTAL  # Arc length is the total line length

if abs(a_catenary) < 1e-6 : # Ans_x was zero or very small
    print("Warning: 'a_catenary' is near zero. Catenary is degenerate (line hangs vertically or drapes).")
    print("Node coordinate calculation for catenary shape might be invalid.")
    # A different geometric model (e.g., vertical line + horizontal on seabed) would be needed.
    # The current script will likely fail or give meaningless results for node positions.
    # For now, exiting if 'a' is too small, as catenary formulas divide by 'a'.
    if abs(L_APFP) < 1e-6: # Vertically aligned
        print("Line hangs vertically. Node positions can be calculated by simple vertical segments.")
        # Node calculation for vertical hang:
        nodes = []
        num_internal_nodes = 19
        num_segments = num_internal_nodes + 1
        s_seg_vert = L_TOTAL / num_segments
        z_start_hang = max(zf,za) # Hangs from the higher point
        x_hang = xf # Same x for all
        y_hang = 0.0
        for i_node in range(1, num_internal_nodes + 1):
            s_node_from_top = i_node * s_seg_vert
            node_z_hang = z_start_hang - s_node_from_top # Assuming z is upwards positive or downwards negative consistently
                                                        # With z downwards positive (seabed at large +z), or z upwards (seabed large -z)
                                                        # Given zf=-70, za=-320, z is depth (more negative = deeper)
                                                        # Higher point is zf = -70.
            node_z_val = zf - s_node_from_top # Nodes go down from Fairlead
            if node_z_val < za and L_TOTAL > h_span : # Passed anchor depth and line has slack
                # This implies part on seabed - more complex.
                # For now, just simple vertical distribution if Ans_x = 0
                pass # This simplified vertical part may not be fully correct for all H=0 scenarios
            nodes.append({"id": i_node, "x": x_hang, "y": y_hang, "z": node_z_val, "s_from_fp": s_node_from_top})
        internal_nodes_coords_final = nodes

    else: # L_APFP > 0, H=0. Complex draped shape.
        print("Ans_x=0 with L_APFP > 0. Catenary formulas invalid. Draped shape not calculated by this script.")
        internal_nodes_coords_final = [] # Cannot calculate with this script.
else:
    # General case: a_catenary > 0
    # Calculate normalized positions x1_n = (xf - x_gm)/a and x2_n = (xa - x_gm)/a
    # sum_norm_x_half = (x1_n + x2_n) / 2
    # diff_norm_x_half = (x2_n - x1_n) / 2 = Xs_rel / (2 * a_catenary)
    # Zs_rel / S_line = tanh(sum_norm_x_half) if sinh(diff_norm_x_half) is used as denominator for Zs_rel/a
    # Zs_rel / S_line = tanh((x1_n+x2_n)/2) (This is a property that relates endpoints and length)
    
    # Robust calculation for x1_n, x2_n, x_gm, z_offset
    # Term for denominator: 2 * sinh(Xs_rel / (2 * a_catenary))
    denominator_term_sum_x = 2 * math.sinh(Xs_rel / (2 * a_catenary))
    if abs(denominator_term_sum_x) < 1e-9: # Avoid division by zero if Xs_rel is small or 'a' is very large
        sum_norm_x_half = Zs_rel / Xs_rel if abs(Xs_rel) > 1e-6 else 0.0 # Approx for nearly straight line
    else:
        sum_norm_x_half = myasinh( (Zs_rel / a_catenary) / denominator_term_sum_x )

    diff_norm_x_half = Xs_rel / (2 * a_catenary)
    
    x1_n = sum_norm_x_half - diff_norm_x_half  # (xf - x_gm)/a
    x2_n = sum_norm_x_half + diff_norm_x_half  # (xa - x_gm)/a

    x_gm = xf - a_catenary * x1_n
    z_offset = zf - a_catenary * math.cosh(x1_n)

    print(f"Global catenary: a={a_catenary:.3f}, x_gm={x_gm:.3f}, z_offset={z_offset:.3f}")
    # print(f"  (xf-xgm)/a = {x1_n:.4f}, (xa-xgm)/a = {x2_n:.4f}")
    # Sanity check with anchor point
    # z_a_recalc = a_catenary * math.cosh(x2_n) + z_offset
    # S_recalc = a_catenary * (math.sinh(x2_n) - math.sinh(x1_n))
    # print(f"  Anchor z check: actual={za:.3f}, recalc={z_a_recalc:.3f}")
    # print(f"  Line length S check: actual={S_line:.3f}, recalc={S_recalc:.3f}")


    # --- Calculate internal node coordinates ---
    num_internal_nodes = 19
    num_segments = num_internal_nodes + 1
    s_segment_len = S_line / num_segments # Arc length of each segment

    internal_nodes_coords_final = []
    # sinh_x1_n = math.sinh(x1_n) # sinh((xf - x_gm)/a)

    for i_node in range(1, num_internal_nodes + 1):
        s_k_arc_from_fp = i_node * s_segment_len # Arc length from Fairlead to current node k
        
        # Solve for xk_n = (node_x_coord - x_gm) / a_catenary
        # s_k_arc_from_fp = a_catenary * (math.sinh(xk_n) - math.sinh(x1_n))
        arg_for_asinh = s_k_arc_from_fp / a_catenary + math.sinh(x1_n)
        xk_n = myasinh(arg_for_asinh)
        
        node_x_coord = x_gm + a_catenary * xk_n
        node_y_coord = 0.0 # As per problem, y=0
        node_z_coord = a_catenary * math.cosh(xk_n) + z_offset
        
        internal_nodes_coords_final.append({
            "id": i_node,
            "x": node_x_coord,
            "y": node_y_coord,
            "z": node_z_coord,
            "s_from_fp": s_k_arc_from_fp
        })

# --- Output Results ---
print("\n--- Calculated 19 Internal Node Coordinates ---")
if not internal_nodes_coords_final and ans_x_final == 0 and abs(L_APFP)>1e-6:
    print("Node coordinates not calculated due to Ans_x = 0 and L_APFP > 0 (draped shape).")
elif not internal_nodes_coords_final and ans_x_final !=0 :
     print("Node coordinates calculation failed or was skipped.")
else:
    print("Fairlead Point: X={:.3f}, Y={:.3f}, Z={:.3f}".format(FP_COORDS['x'], FP_COORDS['y'], FP_COORDS['z']))
    print("----------------------------------------------------------------------")
    print("ID |   X (m)   |   Y (m)   |   Z (m)   | Arc Length from FP (m)")
    print("---|-----------|-----------|-----------|------------------------")
    for node_data in internal_nodes_coords_final:
        print(f"{node_data['id']:<2} | {node_data['x']:9.3f} | {node_data['y']:9.3f} | {node_data['z']:9.3f} | {node_data['s_from_fp']:11.3f}")
    print("----------------------------------------------------------------------")
    print("Anchor Point:   X={:.3f}, Y={:.3f}, Z={:.3f}".format(AP_COORDS['x'], AP_COORDS['y'], AP_COORDS['z']))
