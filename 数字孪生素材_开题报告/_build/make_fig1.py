# -*- coding: utf-8 -*-
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib; matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import _fonts

fig, ax = plt.subplots(figsize=(13.6, 8.0))
ax.set_xlim(0, 100); ax.set_ylim(0, 78); ax.axis('off')
fig.canvas.draw()

texts = {}
def box(cx, cy, w, h, lines, name, fc='white', ec='#31557a', lw=1.5, fs=9.6, tc='#14181c', bold=False):
    x, y = cx-w/2, cy-h/2
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.02,rounding_size=1.2', fc=fc, ec=ec, lw=lw, zorder=3))
    if isinstance(lines, str): lines=[lines]
    texts[name] = (ax.text(cx, cy, '\n'.join(lines), ha='center', va='center', fontsize=fs, color=tc,
                  zorder=5, linespacing=1.35, fontweight=('bold' if bold else 'normal')), cx, cy, w, h)
def arrow(x1, y1, x2, y2, color='#33475b', lw=1.6, style='-|>', ls='-', ms=13, z=2):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style, mutation_scale=ms, lw=lw,
                 color=color, linestyle=ls, zorder=z, shrinkA=2, shrinkB=2))
def arc(x1, y1, x2, y2, rad, color='#b0301f', lw=1.7, style='-|>', ls='--'):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), connectionstyle=f'arc3,rad={rad}',
                 arrowstyle=style, mutation_scale=12, lw=lw, color=color, linestyle=ls, zorder=2))

ax.text(50, 75.5, '面向毫米波雷达无接触血压估计的“胸廓—心血管运动”数字孪生框架',
        ha='center', va='center', fontsize=13.5, fontweight='bold', color='#10233a')

# ========== 物理域（左，橙） ==========
box(12.5, 50, 21, 40, '', name='LBOX')
ax.text(12.5, 67.5, '真实物理域\n(被测对象＋实测传感器)', ha='center', va='center', fontsize=10, color='#8a5a18', fontweight='bold', linespacing=1.4)
tA = box(12.5, 60, 18, 6.4, ['真实人体生理活动','胸壁·心搏·脉搏波·呼吸微动'], 'A', fc='#fae6c4', ec='#a8791f', fs=8.8)
tB = box(12.5, 51.5, 18, 5.8, ['毫米波雷达实测','77 GHz FMCW 回波'], 'B', fc='#f7d9ad', ec='#a8791f', fs=8.8)
tC = box(12.5, 42, 18, 7.0, ['实测信号库＋金标准血压','(SBP/DBP 参考真值)'], 'C', fc='#fae6c4', ec='#a8791f', fs=8.8)
arrow(12.5, 56.7, 12.5, 54.5); arrow(12.5, 48.5, 12.5, 45.6)

# ========== 数字孪生平台（右，蓝容器） ==========
box(62, 41, 68, 52, '', name='PBOX')
ax.text(62, 69.8, '数字孪生仿真平台（虚拟传感与机理化数据生成）', ha='center', va='center',
        fontsize=11.5, color='#1d4e79', fontweight='bold')

tP1 = box(44, 62, 30, 6.0, ['生理参数与个体先验','体型·呼吸节律·心率·血管弹性'], 'P1', fc='#dcebfa', ec='#31557a', fs=8.8)
tP2 = box(84, 62, 22, 6.0, ['测量场景参数','位置·频率·带宽·SNR'], 'P2', fc='#dcebfa', ec='#31557a', fs=8.8)
tM1 = box(32, 51.5, 16.5, 8.6, ['① 几何—解剖模型','胸腔与心脏三维表面','(STL 几何)'], 'M1', fc='#e7edf3', ec='#34506e', fs=8.6)
tM2 = box(56, 51.5, 22, 8.6, ['② 生理运动机理模型','形变包络 A(x)×调制 M(t)','呼吸·心搏·脉动传播'], 'M2', fc='#e0f0e7', ec='#2e6d54', fs=8.4)
tM3 = box(84, 51.5, 16.5, 8.6, ['③ 电磁前向仿真','逐点 FMCW 回波','时延/相位·Lambertian'], 'M3', fc='#ebe3f5', ec='#6b4f9e', fs=8.6)
tO = box(60, 39.5, 46, 8.2, ['④ 虚拟基带差拍信号 Mix ＋ 真值标签','位移场·呼吸/心率·PTT·SBP/DBP → 带标注仿真数据集'], 'O', fc='#e3f1e3', ec='#2f7a3f', fs=8.8)
tS = box(56, 27, 52, 9.4, ['⑤ 模型训练 · 受控评测 · 机理可解释验证','信号分离｜特征提取｜波形重建｜生理特征回归', '跨个体泛化评估 · 呼吸谐波/体位/SNR 消融实验'], 'S', fc='#fff1d6', ec='#c2762a', fs=8.2)
tR = box(89.5, 27, 12.5, 9.4, ['无接触连续血压','与生命体征估计','(SBP/DBP·HR·RR)'], 'R', fc='#fbe3da', ec='#c0562a', fs=8.2)

# ---- 平台内正向流程 ----
arrow(39, 59, 32, 56.1)          # P1 -> M1
arrow(50, 59, 52.5, 56.1)        # P1 -> M2
arrow(84, 59, 84, 56.1)          # P2 -> M3
arrow(40.3, 51.5, 45.5, 51.5)    # M1 -> M2
arrow(67, 51.5, 76, 51.5)        # M2 -> M3
arrow(84, 47.2, 66, 44.1)        # M3 -> O
arrow(60, 35.4, 60, 32.0)        # O -> S
arrow(82.5, 27, 83.6, 27)        # S -> R

# ---- 实测数据流入训练/评测（橙实线） ----
arrow(21.5, 42, 31.5, 29.5, color='#a8791f', lw=1.6)

# ---- 虚实反馈：残差→参数标定/孪生更新（红虚线，右侧大回环） ----
arc(94.5, 24.5, 90.5, 60.5, rad=-0.22)
ax.text(74.5, 71.5, '残差反馈→参数标定·孪生更新·个性化校准', fontsize=8.0, color='#b0301f', ha='center')

ax.text(50, 8.6, '虚实闭环：以实测—仿真对比残差对几何/运动/电磁及个体参数进行标定与更新，使孪生模型与真实生理状态逐步对齐（虚线）。',
        ha='center', va='center', fontsize=9, color='#b0301f')

# ===== 自检 1：文本溢出 =====
renderer = fig.canvas.get_renderer()
def box_px(cx,cy,w,h):
    x,y=cx-w/2,cy-h/2
    (x0,y0)=ax.transData.transform((x,y)); (x1,y1)=ax.transData.transform((x+w,y+h))
    return min(x0,x1),min(y0,y1),max(x0,x1),max(y0,y1)
print('== overflow check ==')
bad=False
for nm,(t,cx,cy,w,h) in texts.items():
    bb=t.get_window_extent(renderer=renderer)
    bx0,by0,bx1,by1=box_px(cx,cy,w,h)
    wr=bb.width/(bx1-bx0); hr=bb.height/(by1-by0)
    if wr>1 or hr>1: bad=True
    print(f'{nm:5s} wr={wr:4.2f} hr={hr:4.2f} {"OVER!" if (wr>1 or hr>1) else ""}')

# ===== 自检 2：框间重叠 =====
boxes = {'A':(12.5,60,18,6.4),'B':(12.5,51.5,18,5.8),'C':(12.5,42,18,7.0),
 'P1':(44,62,30,6),'P2':(84,62,22,6),'M1':(32,51.5,16.5,8.6),'M2':(56,51.5,22,8.6),
 'M3':(84,51.5,16.5,8.6),'O':(60,39.5,46,8.2),'S':(56,27,52,9.4),'R':(89.5,27,12.5,9.4)}
print('== box overlap check ==')
keys=list(boxes); ov=False
for i in range(len(keys)):
    for j in range(i+1,len(keys)):
        a=keys[i]; b=keys[j]
        (cx1,cy1,w1,h1)=boxes[a]; (cx2,cy2,w2,h2)=boxes[b]
        ox=max(0,min(cx1+w1/2,cx2+w2/2)-max(cx1-w1/2,cx2-w2/2))
        oy=max(0,min(cy1+h1/2,cy2+h2/2)-max(cy1-h1/2,cy2-h2/2))
        if ox>0.05 and oy>0.05:
            print(f'OVERLAP {a} & {b}: ox={ox:.2f} oy={oy:.2f}'); ov=True
print('overlap none' if not ov else 'FIX overlaps')
print('overflow none' if not bad else 'FIX overflow')

plt.subplots_adjust(left=0.005, right=0.995, top=0.995, bottom=0.005)
out='数字孪生素材_开题报告/图/fig1_数字孪生总体架构图.png'
fig.savefig(out, dpi=300); fig.savefig(out.replace('.png','.svg'))
print('saved')
