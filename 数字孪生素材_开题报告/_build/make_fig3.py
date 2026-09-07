# -*- coding: utf-8 -*-
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib; matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.patches import FancyArrowPatch
import _fonts

# ---- 复现 wanzheng.m 参数 ----
T_in, T_plateau, T_out = 0.8, 0.4, 0.8
T_b = T_in + T_plateau + T_out + T_plateau   # 2.4 s（呼气末静息段并入周期）
T_h, tau_h, sigma_h = 0.8, 0.8/3, 0.8/5

def M_breath(t):
    tm = np.mod(t, T_b)
    out = np.empty_like(t)
    i1 = tm < T_in
    i2 = (tm >= T_in) & (tm < T_in + T_plateau)
    i3 = (tm >= T_in + T_plateau) & (tm < T_in + T_plateau + T_out)
    i4 = ~(i1|i2|i3)
    out[i1] = np.sin((np.pi/2)*(tm[i1]/T_in))**2
    out[i2] = 1.0
    x = (tm[i3]-T_in-T_plateau)/T_out
    out[i3] = np.sin(np.pi/2 + (np.pi/2)*x)**2
    out[i4] = 0.0
    return out

def M_heart(t):
    tm = np.mod(t, T_h)
    return np.exp(-((tm - tau_h)/sigma_h)**2)

t = np.arange(0, 4.8, 1e-4)
mb, mh = M_breath(t), M_heart(t)
Ar, Ah = 4.0, 0.35          # 呼吸、心搏体表位移幅值示意（mm）
d_breath = Ar*mb
d_heart  = Ah*mh
d_tot    = d_breath + d_heart

fig = plt.figure(figsize=(7.6, 9.2))
gs = fig.add_gridspec(3, 1, height_ratios=[1,1,1.35], hspace=0.42, left=0.10, right=0.97, top=0.965, bottom=0.06)

axa = fig.add_subplot(gs[0])
axa.plot(t, mb, color='#1f77b4', lw=1.8)
axa.set_ylim(-0.08, 1.12); axa.set_xlim(0, 4.8)
axa.set_yticks([0,0.5,1.0]); axa.set_ylabel('呼吸调制幅度 M$_r$(t)')
axa.set_title('(a) 呼吸驱动：吸气—屏气—呼气—静息 四段式（周期 2.4 s）', fontsize=10)
for ph in range(2):
    s = ph*T_b
    for x,c,txt in [(s+0.8,'#d62728','吸气末平台'),(s+2.0,'#2ca02c','呼气末静息')]:
        axa.axvline(x, color=c, lw=0.8, ls=':')
    axa.annotate('吸气',  xy=(s+0.35,0.9), fontsize=8.5, ha='center', color='#1f77b4')
    axa.annotate('呼气',  xy=(s+1.6,0.9), fontsize=8.5, ha='center', color='#9467bd')
axa.grid(alpha=.25, lw=.5)

axb = fig.add_subplot(gs[1])
axb.plot(t, mh, color='#e3772c', lw=1.6)
axb.set_ylim(-0.08, 1.15); axb.set_xlim(0, 4.8)
axb.set_yticks([0,0.5,1.0]); axb.set_ylabel('心搏调制幅度 M$_h$(t)')
axb.set_title('(b) 心搏驱动：高斯脉冲（周期 0.8 s ≈ 75 bpm，τ=T/3，σ=T/5）', fontsize=10)
for k in range(6):
    c = tau_h + k*T_h
    axb.axvline(c, color='#e3772c', lw=0.6, ls=':')
    axb.annotate(f'{c:.2f}s', xy=(c, 1.02), fontsize=6.8, ha='center', color='#8a4a10')
axb.grid(alpha=.25, lw=.5)

axc = fig.add_subplot(gs[2])
axc.plot(t, d_tot, color='#0a3d62', lw=1.4, label='复合体表位移 d(t)=d$_r$+d$_h$')
axc.plot(t, d_breath, color='#9aa5ad', lw=1.0, ls='--', label='呼吸分量 d$_r$(t)（4 mm 示意）')
axc.set_xlim(0, 4.8); axc.set_xlabel('时间 / s'); axc.set_ylabel('体表微动位移 / mm')
axc.set_title('(c) 体表复合微动与心搏分量（呼吸数 mm 级、心搏亚 mm 级，图中相对示意）', fontsize=10)
axc.legend(loc='upper left', fontsize=7.6, frameon=False)
axc.grid(alpha=.25, lw=.5)

# 放大窗口：呼气末静息附近（呼吸分量接近0，心搏脉冲可辨）
x0,x1,y0,y1 = 1.75, 2.1, -0.15, 0.75
axc.axvspan(x0,x1, color='#e3772c', alpha=.10)
axc.annotate('静息段局部放大', xy=((x0+x1)/2, 0.6), xytext=(2.6, 2.2), fontsize=8,
             arrowprops=dict(arrowstyle='-|>', color='#8a4a10', lw=1.0), color='#8a4a10', ha='center')
axz = axc.inset_axes([0.60, 0.50, 0.30, 0.36])
axz.plot(t, d_tot, color='#0a3d62', lw=1.3)
axz.plot(t, d_breath, color='#9aa5ad', lw=.8, ls='--')
axz.set_xlim(x0,x1); axz.set_ylim(y0,y1)
axz.set_xticks([1.8,1.9,2.0]); axz.set_yticks([0,0.4])
axz.set_title('1.75–2.1 s', fontsize=6.5)
axz.grid(alpha=.3, lw=.4)
for c in (1.8667,):
    axz.axvline(c, color='#e3772c', lw=.7, ls=':')

fig.savefig('数字孪生素材_开题报告/图/fig3_呼吸心跳运动规律曲线.png', dpi=300)
fig.savefig('数字孪生素材_开题报告/图/fig3_呼吸心跳运动规律曲线.svg')
print('saved fig3')
