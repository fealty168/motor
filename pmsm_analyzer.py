import streamlit as st
import numpy as np
import pandas as pd
import math

# ---------- 页面设置 ----------
st.set_page_config(page_title="永磁同步电机分析智能体", layout="wide")
st.title("🔧 永磁同步电机 (PMSM) 分析智能体")
st.markdown("根据输入的电机几何、绕组、磁钢参数，自动计算电磁性能与控制器参数。")

# ---------- 侧边栏输入 ----------
st.sidebar.header("📥 输入参数")

with st.sidebar.expander("⚡ 基本电气参数", expanded=True):
    col1, col2 = st.columns(2)
    with col1:
        P = st.number_input("额定功率 (kW)", min_value=0.1, max_value=1000.0, value=75.0)
        V_ll = st.number_input("额定线电压 (V)", min_value=50.0, max_value=1000.0, value=380.0)
        n = st.number_input("额定转速 (rpm)", min_value=100.0, max_value=20000.0, value=3000.0)
    with col2:
        f = st.number_input("额定频率 (Hz)", min_value=10.0, max_value=1000.0, value=200.0)
        I_rated = st.number_input("额定电流 (A)", min_value=1.0, max_value=2000.0, value=122.0)
        poles = st.number_input("极数", min_value=2, max_value=48, value=8, step=2)

with st.sidebar.expander("🧲 磁钢参数", expanded=True):
    col1, col2 = st.columns(2)
    with col1:
        mag_thick = st.number_input("磁钢厚度 (mm)", min_value=1.0, max_value=20.0, value=6.0) / 1000
        mag_width = st.number_input("单块磁钢宽度 (mm)", min_value=5.0, max_value=100.0, value=38.0) / 1000
        mag_length = st.number_input("磁钢长度 (mm)", min_value=10.0, max_value=500.0, value=120.0) / 1000
        per_pole_pieces = st.number_input("每极磁钢片数", min_value=1, max_value=4, value=2)
    with col2:
        Br = st.number_input("剩磁 Br (T)", min_value=0.5, max_value=1.5, value=1.24)
        mu_r = st.number_input("相对磁导率 μr", min_value=1.0, max_value=1.5, value=1.05)
        bridge = st.number_input("磁桥厚度 (mm)", min_value=0.5, max_value=5.0, value=2.5) / 1000

with st.sidebar.expander("📐 定子几何", expanded=True):
    col1, col2 = st.columns(2)
    with col1:
        D_out = st.number_input("定子外径 (mm)", min_value=50.0, max_value=1000.0, value=368.0) / 1000
        D_in = st.number_input("定子内径 (mm)", min_value=30.0, max_value=900.0, value=260.0) / 1000
        L = st.number_input("铁芯长度 (mm)", min_value=20.0, max_value=500.0, value=120.0) / 1000
        airgap = st.number_input("气隙 (mm)", min_value=0.2, max_value=5.0, value=1.1) / 1000
        slots = st.number_input("槽数", min_value=6, max_value=96, value=48)
    with col2:
        Bs0 = st.number_input("槽口宽度 (mm)", min_value=1.0, max_value=10.0, value=3.2) / 1000
        Hs0 = st.number_input("槽口高度 (mm)", min_value=0.2, max_value=5.0, value=1.0) / 1000
        Hs1 = st.number_input("槽肩高度 (mm)", min_value=0.2, max_value=10.0, value=1.26) / 1000
        Hs2 = st.number_input("槽深 (mm)", min_value=5.0, max_value=50.0, value=25.64) / 1000
        Bs2_radius = st.number_input("槽底圆角半径 (mm)", min_value=1.0, max_value=20.0, value=4.91) / 1000
        tooth_shoulder = st.number_input("槽肩宽度 (mm)", min_value=2.0, max_value=20.0, value=6.5) / 1000

with st.sidebar.expander("🧶 绕组参数", expanded=True):
    col1, col2 = st.columns(2)
    with col1:
        turns_per_coil = st.number_input("每线圈匝数", min_value=1, max_value=100, value=9)
        layers = st.number_input("层数", min_value=1, max_value=2, value=2)
        parallel_branches = st.number_input("并联支路数", min_value=1, max_value=8, value=4)
        coil_span = st.number_input("线圈跨距 (槽数)", min_value=1, max_value=slots//2, value=5)
    with col2:
        wire_dia = st.number_input("裸线径 (mm)", min_value=0.1, max_value=5.0, value=0.95) / 1000
        ins_dia = st.number_input("含漆层直径 (mm)", min_value=0.1, max_value=5.0, value=1.01) / 1000
        parallel_wires = st.number_input("并绕根数", min_value=1, max_value=100, value=10)
        end_length = st.number_input("端部平均半匝长 (mm)", min_value=50.0, max_value=1000.0, value=341.45) / 1000

with st.sidebar.expander("📊 材料与系数", expanded=True):
    col1, col2 = st.columns(2)
    with col1:
        Kc = st.number_input("卡特系数", min_value=1.0, max_value=2.0, value=1.1)
        Ksat = st.number_input("饱和系数", min_value=1.0, max_value=2.0, value=1.1)
        sigma = st.number_input("漏磁系数", min_value=1.0, max_value=1.5, value=1.1)
        alpha_i = st.number_input("计算极弧系数", min_value=0.5, max_value=1.0, value=0.75)
    with col2:
        rho = st.number_input("电阻率 @20°C (Ω·mm²/m)", min_value=0.01, max_value=0.05, value=0.0217)
        temp_coeff = st.number_input("温度系数 (/°C)", min_value=0.0, max_value=0.01, value=0.0039)
        temp_rise = st.number_input("温升 (°C)", min_value=0, max_value=150, value=55)

with st.sidebar.expander("📏 实测数据覆盖（可选）", expanded=False):
    use_measured = st.checkbox("使用实测电感覆盖计算值")
    if use_measured:
        Ld_meas = st.number_input("实测 Ld (mH)", min_value=0.01, max_value=100.0, value=0.5205) / 1000
        Lq_meas = st.number_input("实测 Lq (mH)", min_value=0.01, max_value=100.0, value=1.224) / 1000

# ---------- 核心计算函数 ----------
def calculate(p):
    # 基本常数
    p_pole_pairs = p['poles'] / 2
    # 极距
    tau = math.pi * p['D_in'] / p['poles']  # 米
    # 转子外径
    D_rotor = p['D_in'] - 2 * p['airgap']
    # 每槽导体数
    conductors_per_slot = p['turns_per_coil'] * p['layers']
    # 总导体数
    total_conductors = conductors_per_slot * p['slots']
    # 每相串联匝数
    Nph = total_conductors / (3 * p['parallel_branches'])
    # 绕组系数 (简化计算：短距系数 * 分布系数)
    # 分布系数: q=slots/3/poles, 对于整数槽 q=2 (假设48槽8极)
    q = p['slots'] / (3 * p['poles'])
    if abs(q - 2) < 0.1:
        Kd = 0.966  # 常见值
    else:
        Kd = math.sin(math.pi/6) / (q * math.sin(math.pi/6/q))  # 通用公式，此处简化
    # 短距系数
    full_pitch = p['slots'] / p['poles']  # 整距槽数
    beta = p['coil_span'] / full_pitch
    Kp = math.sin(beta * math.pi/2)
    Kw = Kd * Kp

    # ---------- 磁路估算 ----------
    # d轴等效气隙 (考虑磁钢)
    g_equiv_d = p['airgap'] * p['Kc'] * p['Ksat'] + p['mag_thick'] / p['mu_r']
    # 气隙磁密基波幅值 (磁路法)
    Bg1 = p['Br'] * (p['mag_thick'] / p['mu_r']) / (p['mag_thick'] / p['mu_r'] + p['airgap'] * p['Kc'] * p['Ksat']) * p['alpha_i'] / p['sigma']
    # 每极磁通
    Phi = (2/math.pi) * Bg1 * tau * p['L']
    # 空载相反电势 (有效值)
    E_ph = 4.44 * p['f'] * Nph * Phi * Kw
    # 线反电势
    E_line = E_ph * math.sqrt(3)

    # ---------- 永磁磁链 ----------
    lambda_f = E_ph * math.sqrt(2) / (2 * math.pi * p['f'])

    # ---------- 电感估算 ----------
    # 槽漏感等先忽略，采用磁路法估算主电感
    # (Nph*Kw)^2
    NK_sq = (Nph * Kw) ** 2
    # d轴主电感
    Lmd = (3/math.pi) * 4e-7 * math.pi * NK_sq * tau * p['L'] / (p_pole_pairs * g_equiv_d)
    # q轴等效气隙 (无磁钢)
    g_equiv_q = p['airgap'] * p['Kc'] * p['Ksat']
    Lmq = (3/math.pi) * 4e-7 * math.pi * NK_sq * tau * p['L'] / (p_pole_pairs * g_equiv_q)
    # 漏感估算 (取0.15 * Lmq)
    Lsigma = 0.15 * Lmq
    Ld_calc = Lmd + Lsigma
    Lq_calc = Lmq + Lsigma

    # 若用户输入实测电感，则覆盖
    if p['use_measured']:
        Ld = p['Ld_meas']
        Lq = p['Lq_meas']
    else:
        Ld = Ld_calc
        Lq = Lq_calc

    # ---------- 电阻估算 ----------
    # 平均半匝长: 端部长度 + 铁芯长
    avg_half_turn = p['L'] + p['end_length']
    total_cond_length = Nph * 2 * avg_half_turn * p['parallel_branches']  # 每相总长
    # 单根导线截面积
    wire_area = math.pi * (p['wire_dia']/2)**2
    # 总截面积 (每支路)
    total_wire_area = wire_area * p['parallel_wires']
    # 20°C电阻
    R20 = p['rho'] * total_cond_length / total_wire_area
    # 热态电阻
    R_hot = R20 * (1 + p['temp_coeff'] * p['temp_rise'])

    # ---------- 槽满率 ----------
    slot_area = (p['Bs2_radius']*2 + p['tooth_shoulder'])/2 * p['Hs2']  # 简化梯形
    if slot_area < 1e-6:
        slot_area = 200e-6  # 默认
    # 每槽导线总面积 (含绝缘)
    wire_ins_area = math.pi * (p['ins_dia']/2)**2
    copper_area = conductors_per_slot * p['parallel_wires'] * wire_ins_area
    fill_factor = copper_area / slot_area

    # ---------- 电流密度 ----------
    I_phase = p['I_rated']  # 星形线电流等于相电流
    I_per_branch = I_phase / p['parallel_branches']
    J = I_per_branch / total_wire_area  # A/mm²

    # ---------- 性能估算 ----------
    # i_d=0 转矩
    I_peak = I_phase * math.sqrt(2)
    Te_id0 = 1.5 * p_pole_pairs * lambda_f * I_peak
    # 额定转矩
    T_rated = p['P'] * 1000 / (2 * math.pi * p['n'] / 60)

    # 特征电流
    I_ch_peak = lambda_f / Ld
    I_ch_rms = I_ch_peak / math.sqrt(2)

    # 转折速度 (空载反电势等于额定电压)
    n_base = p['n'] * (p['V_ll'] / E_line) if E_line > 0 else p['n']

    # 弱磁估算 (近似)
    # 在额定电压下，输出额定转矩所需弱磁电流简化计算
    # 电压极限方程 (忽略电阻)
    # V_s^2 = (ω λ_f + ω Ld Id)^2 + (ω Lq Iq)^2
    # 此处简化处理，仅提示是否需要弱磁
    V_ph_peak = p['V_ll'] * math.sqrt(2) / math.sqrt(3)
    omega = 2 * math.pi * p['f']
    # 在i_d=0下所需电压
    V_id0 = omega * math.sqrt((lambda_f)**2 + (Lq * I_peak)**2)
    if V_id0 <= V_ph_peak:
        need_flux_weakening = "否 (i_d=0可行)"
        Id_rms_needed = 0.0
    else:
        need_flux_weakening = "是 (需负d轴电流)"
        # 粗略估算所需Id (忽略电阻)
        # 解方程 (λ_f + Ld Id)^2 + (Lq Iq)^2 = (V_ph_peak/ω)^2, 且 Id^2+Iq^2=I_peak^2
        # 近似：令Iq=I_peak, 求所需Id
        V_limit = V_ph_peak / omega
        # (λ_f + Ld Id) ≈ sqrt(V_limit^2 - (Lq I_peak)^2) 或直接近似
        try:
            Id_peak_needed = (math.sqrt(V_limit**2 - (Lq * I_peak)**2) - lambda_f) / Ld
        except:
            Id_peak_needed = -I_peak * 0.5  # 默认
        Id_rms_needed = Id_peak_needed / math.sqrt(2)

    # 凸极率
    saliency = Lq / Ld

    # ---------- 组装结果 ----------
    results = {
        "每相串联匝数": Nph,
        "绕组系数 Kw": Kw,
        "气隙磁密基波 (T)": Bg1,
        "每极磁通 (Wb)": Phi,
        "空载相反电势 (V)": E_ph,
        "空载线反电势 (V)": E_line,
        "永磁磁链 λ_f (Wb)": lambda_f,
        "d轴电感 Ld (mH)": Ld * 1000,
        "q轴电感 Lq (mH)": Lq * 1000,
        "凸极率 Lq/Ld": saliency,
        "相电阻 (热态, Ω)": R_hot,
        "槽满率 (%)": fill_factor * 100,
        "电流密度 (A/mm²)": J,
        "i_d=0 转矩 (N·m)": Te_id0,
        "额定转矩 (N·m)": T_rated,
        "特征电流 (有效值, A)": I_ch_rms,
        "转折速度 (rpm)": n_base,
        "需要弱磁": need_flux_weakening,
        "估算d轴弱磁电流 (有效值, A)": Id_rms_needed,
    }
    return results

# ---------- 执行计算 ----------
params = {
    'P': P, 'V_ll': V_ll, 'n': n, 'f': f, 'I_rated': I_rated, 'poles': poles,
    'mag_thick': mag_thick, 'mag_width': mag_width, 'mag_length': mag_length,
    'per_pole_pieces': per_pole_pieces, 'Br': Br, 'mu_r': mu_r, 'bridge': bridge,
    'D_out': D_out, 'D_in': D_in, 'L': L, 'airgap': airgap, 'slots': slots,
    'Bs0': Bs0, 'Hs0': Hs0, 'Hs1': Hs1, 'Hs2': Hs2, 'Bs2_radius': Bs2_radius,
    'tooth_shoulder': tooth_shoulder,
    'turns_per_coil': turns_per_coil, 'layers': layers, 'parallel_branches': parallel_branches,
    'coil_span': coil_span, 'wire_dia': wire_dia, 'ins_dia': ins_dia,
    'parallel_wires': parallel_wires, 'end_length': end_length,
    'Kc': Kc, 'Ksat': Ksat, 'sigma': sigma, 'alpha_i': alpha_i,
    'rho': rho, 'temp_coeff': temp_coeff, 'temp_rise': temp_rise,
    'use_measured': use_measured,
}
if use_measured:
    params['Ld_meas'] = Ld_meas
    params['Lq_meas'] = Lq_meas

results = calculate(params)

# ---------- 显示结果 ----------
st.header("📊 计算结果")

col1, col2, col3 = st.columns(3)
with col1:
    st.metric("每相串联匝数", f"{results['每相串联匝数']:.1f}")
    st.metric("绕组系数", f"{results['绕组系数 Kw']:.3f}")
    st.metric("气隙磁密基波", f"{results['气隙磁密基波 (T)']:.3f} T")
with col2:
    st.metric("空载相反电势", f"{results['空载相反电势 (V)']:.1f} V")
    st.metric("空载线反电势", f"{results['空载线反电势 (V)']:.1f} V")
    st.metric("永磁磁链", f"{results['永磁磁链 λ_f (Wb)']:.4f} Wb")
with col3:
    st.metric("d轴电感", f"{results['d轴电感 Ld (mH)']:.4f} mH")
    st.metric("q轴电感", f"{results['q轴电感 Lq (mH)']:.4f} mH")
    st.metric("凸极率", f"{results['凸极率 Lq/Ld']:.2f}")

st.subheader("🔋 性能指标")
col1, col2, col3 = st.columns(3)
with col1:
    st.metric("电流密度", f"{results['电流密度 (A/mm²)']:.2f} A/mm²")
    st.metric("槽满率", f"{results['槽满率 (%)']:.1f} %")
    st.metric("相电阻 (热态)", f"{results['相电阻 (热态, Ω)']:.4f} Ω")
with col2:
    st.metric("i_d=0 转矩", f"{results['i_d=0 转矩 (N·m)']:.1f} N·m")
    st.metric("额定转矩", f"{results['额定转矩 (N·m)']:.1f} N·m")
    st.metric("特征电流", f"{results['特征电流 (有效值, A)']:.1f} A")
with col3:
    st.metric("转折速度", f"{results['转折速度 (rpm)']:.0f} rpm")
    st.metric("需要弱磁", results['需要弱磁'])
    st.metric("估算弱磁电流", f"{results['估算d轴弱磁电流 (有效值, A)']:.1f} A")

st.subheader("⚙️ 控制器参数表")
ctrl_df = pd.DataFrame({
    "参数": ["永磁磁链 λ_f (Wb)", "d轴电感 Ld (mH)", "q轴电感 Lq (mH)", "相电阻 (Ω, 热态)",
             "反电势常数 Ke (V_rms/rpm)", "转矩常数 Kt (N·m/A_rms, i_d=0)", "特征电流 (A_rms)"],
    "数值": [f"{results['永磁磁链 λ_f (Wb)']:.4f}",
             f"{results['d轴电感 Ld (mH)']:.4f}",
             f"{results['q轴电感 Lq (mH)']:.4f}",
             f"{results['相电阻 (热态, Ω)']:.4f}",
             f"{results['空载相反电势 (V)']/n:.4f}",
             f"{results['i_d=0 转矩 (N·m)'] / I_rated:.2f}",
             f"{results['特征电流 (有效值, A)']:.1f}"]
})
st.table(ctrl_df)

# ---------- 设计评估 ----------
st.subheader("📌 设计评估")
warnings = []
if results['槽满率 (%)'] > 80:
    warnings.append("⚠️ 槽满率过高 (>80%)，建议减小线径或减少并绕根数。")
elif results['槽满率 (%)'] > 75:
    warnings.append("🟡 槽满率略高 (75-80%)，注意下线工艺。")
else:
    warnings.append("✅ 槽满率合适。")

if results['电流密度 (A/mm²)'] > 8:
    warnings.append("⚠️ 电流密度过高 (>8 A/mm²)，需加强冷却或增大导线截面。")
elif results['电流密度 (A/mm²)'] > 6:
    warnings.append("🟡 电流密度略高 (6-8 A/mm²)，建议改善散热。")
else:
    warnings.append("✅ 电流密度合适。")

if results['空载线反电势 (V)'] > 1.2 * V_ll:
    warnings.append(f"⚠️ 空载反电势 ({results['空载线反电势 (V)']:.0f}V) 超过额定电压1.2倍，弱磁需求大。")
elif results['空载线反电势 (V)'] < 0.8 * V_ll:
    warnings.append(f"🟡 空载反电势 ({results['空载线反电势 (V)']:.0f}V) 偏低，可能影响转矩密度。")
else:
    warnings.append("✅ 反电势匹配良好。")

if results['需要弱磁'] == "是 (需负d轴电流)":
    warnings.append(f"🟡 额定点需弱磁，估算 d 轴电流 {results['估算d轴弱磁电流 (有效值, A)']:.1f} A。")

for w in warnings:
    st.markdown(w)

st.markdown("---")
st.caption("注：计算结果基于解析公式，实际设计请用有限元验证。")