# -*- coding: utf-8 -*-
# 把 fig3_呼吸心跳运动规律曲线 的 (a)(b)(c) 三块分别输出为三张独立图
# 样式/参数与 make_fig3.py 完全一致（同样可在 new/ 目录下运行本脚本）
import sys, os
_build = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _build)
import numpy as np
import matplotlib; matplotlib.use('Agg')
from matplotlib import pyplot as plt
import _fonts

newdir = os.path.dirname(_build)                      # _build 的上级 = 数字孪生素材_开题报告
outdir  = os.path.join(newdir, '图')
os.makedirs(outdir, exist_ok=True)

# ---- 复现 wanzheng.m 参数 ----
T_in, T_plateau, T_out = 0.8, 0.4, 0.8
T_b = T_in + T_plateau + T_out + T_plateau   # 2.4 s
T_h, tau_h, sigma_h = 0.8, 0.8/3, 0.8/5

def M_breath(t):
    tm = np.mod(t, T_b)
    out = np.empty_like(t)
    i1 = tm < T_in; i2 = (tm >= T_in) & (tm < T_in + T_plateau)
    i3 = (tm >= T_in + T_plateau) & (tm < T_in + T_plateau + T_out); i4 = ~(i1|i2|i3)
    out[i1] = np.sin((np.pi/2)*(tm[i1]/T_in))**2
    out[i2] = 1.0
    x = (tm[i3]-T_in-T_plateau)/T_out
    out[i3] = np.sin(np.pi/2 + (np.pi/2)*x)**2
    out[i4] = 0.0
    return out

def M_heart(t):
    return np.exp(-((np.mod(t, T_h) - tau_h)/sigma_h)**2)

t = np.arange(0, 4.8, 1e-4)
mb, mh = M_breath(t), M_heart(t)
Ar, Ah = 4.0, 0.35
d_breath = Ar*mb; d_heart = Ah*mh; d_tot = d_breath + d_heart

def base_ax(title):
    fig, ax = plt.subplots(figsize=(7.6, 3.1), dpi=100)
    ax.set_title(title, fontsize=11)
    ax.grid(alpha=.25, lw=.5)
    return fig, ax

# ---- (a) 呼吸四段式 ----
fig, axa = base_ax('呼吸驱动：吸气—屏气—呼气—静息 四段式（周期 2.4 s）')
axa.plot(t, mb, color='#1f77b4', lw=1.8)
axa.set_ylim(-0.08, 1.12); axa.set_xlim(0, 4.8)
axa.set_yticks([0,0.5,1.0]); axa.set_ylabel('呼吸调制幅度 M$_r$(t)'); axa.set_xlabel('时间 / s')
for ph in range(2):
    s = ph*T_b
    for x,c,txt in [(s+0.8,'#d62728','吸气末平台'),(s+2.0,'#2ca02c','呼气末静息')]:
        axa.axvline(x, color=c, lw=0.8, ls=':')
    axa.annotate('吸气', xy=(s+0.35,0.9), fontsize=8.5, ha='center', color='#1f77b4')
    axa.annotate('呼气', xy=(s+1.6,0.9), fontsize=8.5, ha='center', color='#9467bd')
fig.tight_layout()
for ext in ('png','svg'):
    fig.savefig(os.path.join(outdir, f'fig3a_呼吸驱动四段式曲线.{ext}'), dpi=300 if ext=='png' else None)
plt.close(fig)

# ---- (b) 心搏高斯脉冲 ----
fig, axb = base_ax('心搏驱动：高斯脉冲（周期 0.8 s ≈ 75 bpm，τ=T/3，σ=T/5）')
axb.plot(t, mh, color='#e3772c', lw=1.6)
axb.set_ylim(-0.08, 1.15); axb.set_xlim(0, 4.8)
axb.set_yticks([0,0.5,1.0]); axb.set_ylabel('心搏调制幅度 M$_h$(t)'); axb.set_xlabel('时间 / s')
for k in range(6):
    c = tau_h + k*T_h
    axb.axvline(c, color='#e3772c', lw=0.6, ls=':')
    axb.annotate(f'{c:.2f}s', xy=(c, 1.02), fontsize=6.8, ha='center', color='#8a4a10')
fig.tight_layout()
for ext in ('png','svg'):
    fig.savefig(os.path.join(outdir, f'fig3b_心搏驱动高斯脉冲曲线.{ext}'), dpi=300 if ext=='png' else None)
plt.close(fig)

# ---- (c) 体表复合微动（含静息段放大，同原 fig3 面板 c） ----
fig, axc = plt.subplots(figsize=(7.6, 3.9), dpi=100)
axc.plot(t, d_tot, color='#0a3d62', lw=1.4, label='复合体表位移 d(t)=d$_r$+d$_h$')
axc.plot(t, d_breath, color='#9aa5ad', lw=1.0, ls='--', label='呼吸分量 d$_r$(t)（4 mm 示意）')
axc.set_xlim(0, 4.8); axc.set_xlabel('时间 / s'); axc.set_ylabel('体表微动位移 / mm')
axc.set_title('体表复合微动与心搏分量（呼吸数 mm 级、心搏亚 mm 级，图中相对示意）', fontsize=11)
axc.legend(loc='upper left', fontsize=7.6, frameon=False)
axc.grid(alpha=.25, lw=.5)
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
fig.tight_layout()
for ext in ('png','svg'):
    fig.savefig(os.path.join(outdir, f'fig3c_体表复合微动曲线.{ext}'), dpi=300 if ext=='png' else None)
plt.close(fig)

print('saved fig3a/fig3b/fig3c')
