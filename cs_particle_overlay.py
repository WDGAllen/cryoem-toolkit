#!/usr/bin/env python3
"""CryoSPARC 2D-class particle overlay visualisation."""
import argparse, csv, re
from pathlib import Path
import numpy as np
from PIL import Image, ImageEnhance, ImageColor
from scipy.ndimage import rotate, binary_dilation, binary_erosion, label, gaussian_filter
from scipy.fft import fft2, ifft2, fftfreq
import mrcfile

def read_star(path):
    lines = Path(path).read_text().splitlines()
    start = next(i for i,l in enumerate(lines) if l.strip() == 'data_particles')
    loop = next(i for i in range(start, len(lines)) if lines[i].strip() == 'loop_')
    heads=[]; i=loop+1
    while i<len(lines) and lines[i].strip().startswith('_'):
        heads.append(lines[i].split()[0]); i+=1
    rows=[]
    for l in lines[i:]:
        if not l.strip() or l.lstrip().startswith('#') or l.startswith('data_') or l.strip()=='loop_': continue
        vals=l.split()
        if len(vals) >= len(heads): rows.append(dict(zip(heads, vals)))
    return rows

def basename(s): return Path(s).name
def mic_key(s):
    n=basename(s)
    n=re.sub(r'_thumb_@(?:1x|2x)\.png$','',n)
    n=re.sub(r'\.mrc(?:s)?$','',n)
    # CryoSPARC may prepend a new blob UID and preprocessing may append
    # dose-weighting/denoising labels. The stable cross-job identifier is the
    # FoilHole acquisition ID before the processing suffix.
    m=re.search(r'(FoilHole_.+?)(?:_fractions_patch_aligned|$)',n)
    if m: return m.group(1)
    return n

def class_map(path):
    rows=read_star(path); out={}
    for r in rows:
        n=int(r['_rlnImageName'].split('@')[0])
        out[n]=int(r['_rlnImageName'].split('@')[0])
    # templates_sel.star is a selected subset: key is original class ID,
    # stack index is the order in its data_particles loop.
    out={int(r['_rlnImageName'].split('@')[0]): i+1 for i,r in enumerate(rows)}
    return out

def load_stack(path):
    with mrcfile.open(path, permissive=True) as m: return np.array(m.data, dtype=np.float32)

def lowpass(a, pixel, resolution):
    fy=fftfreq(a.shape[0]); fx=fftfreq(a.shape[1]); yy,xx=np.meshgrid(fy,fx,indexing='ij')
    mask=(xx*xx+yy*yy) <= (pixel/resolution)**2
    return np.real(ifft2(fft2(a)*mask)).astype(np.float32)

def autocontrast(a):
    lo,hi=np.percentile(a,[1,99]); return np.clip((a-lo)/(hi-lo+1e-8),0,1)

def load_cs(path): return np.load(path,allow_pickle=True)

COMMON_COLORS = """Maroon #800000\nRed #FF0000\nPink #FFC0CB\nOrange #FFA500\nKhaki #F0E68C\nYellow #FFFF00\nBlue #0000FF\nLightBlue #ADD8E6\nGreen #008000\nLightGreen #90EE90"""

def parse_color(value):
    v=value.strip().lower().replace('_','').replace('-','').replace(' ','')
    if not v.startswith('#'):
        value=v
    try: return ImageColor.getrgb(value)
    except ValueError: raise ValueError(f"Unknown HTML/X11 colour '{value}'")

def main():
    ap=argparse.ArgumentParser(description='CryoSPARC-native ReconSil particle overlay visualiser.',epilog='''\nBrief usage:\n  python3 cs_particle_overlay.py JOB [STAR COLOUR [MODE] ...] [JOB [STAR COLOUR [MODE] ...] ...] --mics MICS_JOB\n\nEach select-2D JOB is followed by zero or more particle-subset STAR files. A STAR may be followed by a colour and optional mode: border or border dashed. Colours may be #RRGGBB hex values or any HTML/X11 colour name (case, spaces and hyphens are ignored). Without a colour, the subset remains grayscale.\n\nFrequently used colours:\n'''+COMMON_COLORS+'\n\nDetailed outputs:\n  With subsets: original, combined repositioned overlay, then one transparent template-only PNG per STAR subset, followed by all coloured templates.\n  Without subsets: repositioned, original, and transparent templates PNGs.\n  Default mode recolours the combined template. border leaves it grayscale and adds a 4-pixel coloured edge; border dashed makes that edge dashed. Transparent subset outputs remain grayscale.''',formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('inputs',nargs='+',help='select-2D jobs and STAR/colour subset pairs; repeat job groups')
    ap.add_argument('--mics',required=True,help='micrograph-thumbnail generation job')
    ap.add_argument('--ori-dimension',type=float,default=5760,help='long dimension of the original full-resolution micrograph')
    ap.add_argument('--Apix',type=float,default=.85,help='original micrograph pixel size in Angstroms')
    ap.add_argument('--lowpassA',type=float,default=None,help='optional low-pass resolution in Angstroms')
    ap.add_argument('--out',default='J6319/reconsil_raw')
    ap.add_argument('--opacity',type=float,default=1.0); ap.add_argument('--no-autocontrast',action='store_true')
    ap.add_argument('--no-invert-classes',action='store_true', help='Do not invert class-average contrast')
    ap.add_argument('--mask-threshold',type=float,default=.04, help='threshold in original template MRC units')
    ap.add_argument('--mask-pad',type=int,default=10, help='binary dilation radius in template pixels')
    ap.add_argument('--border-thickness',type=int,default=4, help='border thickness in output pixels (5 px original × 0.75, rounded)')
    ap.add_argument('--mask-smoothing',type=float,default=1.5, help='Gaussian smoothing sigma for mask edges in output pixels')
    args=ap.parse_args(); out=Path(args.out); out.mkdir(parents=True,exist_ok=True)
    tj=Path(args.mics); thumb_cs_path=next(tj.glob('*_micrograph_thumbnails.cs')); thumbs_cs=load_cs(thumb_cs_path)
    jobs=[]; current=None; subset_specs=[]
    i=0
    while i < len(args.inputs):
        item=args.inputs[i]
        if item.lower().endswith('.star'):
            if current is None: ap.error('subset STAR must follow a select-2D job')
            colour=None; mode='recolour'; consumed=1
            if i+1<len(args.inputs) and not args.inputs[i+1].lower().endswith('.star'):
                try: colour=parse_color(args.inputs[i+1]); consumed=2
                except ValueError:
                    # The next token is another job name: this STAR has no
                    # colour and remains grayscale.
                    colour=None
            if colour is not None and i+consumed<len(args.inputs) and args.inputs[i+consumed].lower()=='border':
                mode='border'; consumed+=1
                if i+consumed<len(args.inputs) and args.inputs[i+consumed].lower()=='dashed': mode='border_dashed'; consumed+=1
            current[1].append((item,colour,mode)); subset_specs.append((item,colour,mode)); i+=consumed; continue
        else:
            current=[Path(item),[]]; jobs.append(current)
        i+=1
    allowed_sets=[]; subset_labels=[]; subset_colors=[]; subset_modes=[]
    for sj,subset_stars in jobs:
        for subset_star,colour,mode in subset_stars:
            allowed=set()
            for sr in read_star(subset_star):
                image_name=sr.get('_rlnImageName','')
                if '@' in image_name:
                    n,path=image_name.split('@',1)
                    allowed.add((Path(path).name,int(n)))
            allowed_sets.append(allowed); subset_labels.append(Path(subset_star).stem.replace('_particles','')); subset_colors.append(colour); subset_modes.append(mode)
    source_data=[]
    for source_index,(sj,subset_stars) in enumerate(jobs):
        particle_path=sj/'particles_selected.cs'; pass_path=sj/f'{sj.name}_passthrough_particles_selected.cs'
        template_path=sj/'templates_selected.cs'; stack_path=sj/'templates_selected.mrc'
        particles=load_cs(particle_path); passthrough=load_cs(pass_path); cs=load_cs(template_path)
        cmap={int(r['blob/idx'])+1:int(r['blob_selected/idx'])+1 for r in cs}
        source_data.append((particles,passthrough,cmap,load_stack(stack_path),len(subset_stars)))
    thumbs={}
    for tr in thumbs_cs:
        tp=Path(tr['micrograph_thumbnail_blob_1x/path'].decode())
        if tp.exists(): thumbs[mic_key(tp.name)]=tp
    rows=[]; subset_offset=0
    for source_index,(particles,passthrough,cmap,stack,nsub) in enumerate(source_data):
        source_allowed=allowed_sets[subset_offset:subset_offset+nsub]
        locations={int(x['uid']):x for x in passthrough}
        for r in particles:
            blob_path=r['blob/path'].decode() if isinstance(r['blob/path'],bytes) else str(r['blob/path'])
            particle_key=(Path(blob_path).name,int(r['blob/idx'])+1)
            subset_indices=[subset_offset+i for i,s in enumerate(source_allowed) if particle_key in s]
            if allowed_sets and not subset_indices: continue
            loc=locations.get(int(r['uid']))
            if loc is None: continue
            shape=loc['location/micrograph_shape']; mic=loc['location/micrograph_path'].decode()
            mic_psize=float(loc['location/micrograph_psize_A']); align_psize=float(r['alignments2D/psize_A'])
            sx=float(r['alignments2D/shift'][0])*align_psize/mic_psize; sy=float(r['alignments2D/shift'][1])*align_psize/mic_psize
            rows.append({'_rlnMicrographName':mic,'_rlnClassNumber':str(int(r['alignments2D/class'])+1),'_rlnAnglePsi':str(np.degrees(float(r['alignments2D/pose']))),'_x_frac':float(loc['location/center_x_frac'])-sx/float(shape[1]),'_y_frac':float(loc['location/center_y_frac'])-sy/float(shape[0]),'_subset_indices':subset_indices,'_source_index':source_index})
        subset_offset+=nsub
    """legacy block removed"""
    '''
    for subset_star in args.subset_stars:
        allowed=set()
        for sr in read_star(subset_star):
            image_name=sr.get('_rlnImageName','')
            if '@' in image_name:
                n,path=image_name.split('@',1)
            allowed.add((Path(path).name,int(n)))
        allowed_sets.append(allowed)
    # rlnClassNumber is the original 1-based class number. CryoSPARC stores
    # the corresponding original 0-based image index in blob/idx and the
    # selected-stack 0-based index in blob_selected/idx.
    cmap={int(r['blob/idx'])+1:int(r['blob_selected/idx'])+1 for r in cs}
    stack=load_stack(stack_path)
    thumbs={}
    for r in thumbs_cs:
        tp=Path(r['micrograph_thumbnail_blob_1x/path'].decode())
        if tp.exists(): thumbs[mic_key(tp.name)]=tp
    # Join selected-2D alignment data to passthrough exposure/location data.
    locations={int(r['uid']):(r['location/micrograph_path'].decode(),r) for r in passthrough}
    rows=[]
    for r in particles:
        blob_path=r['blob/path'].decode() if isinstance(r['blob/path'],bytes) else str(r['blob/path'])
        particle_key=(Path(blob_path).name,int(r['blob/idx'])+1)
        subset_indices=[i for i,s in enumerate(allowed_sets) if particle_key in s]
        if allowed_sets and not subset_indices:
            continue
        if not allowed_sets:
            subset_indices=[]
        loc=locations.get(int(r['uid']))
        if loc is None: continue
        mic,er=loc; shape=er['location/micrograph_shape']
        mic_psize=float(er['location/micrograph_psize_A'])
        align_psize=float(r['alignments2D/psize_A'])
        sx=float(r['alignments2D/shift'][0])*align_psize/mic_psize
        sy=float(r['alignments2D/shift'][1])*align_psize/mic_psize
        rows.append({'_rlnMicrographName':mic,'_rlnClassNumber':str(int(r['alignments2D/class'])+1),
                     '_rlnAnglePsi':str(np.degrees(float(r['alignments2D/pose']))),
                     '_x_frac':float(er['location/center_x_frac'])-sx/float(shape[1]),
                     '_y_frac':float(er['location/center_y_frac'])-sy/float(shape[0]),
                     '_subset_indices':subset_indices})
    '''
    grouped={}
    for r in rows:
        k=mic_key(r['_rlnMicrographName'])
        if k in thumbs: grouped.setdefault(k,[]).append(r)
    report=[]
    for k,prs in grouped.items():
        p=thumbs[k]; im=np.asarray(Image.open(p).convert('F'),dtype=np.float32)
        thumb_long=max(im.shape); scale=thumb_long/args.ori_dimension
        effective_pixel=args.Apix*scale
        if args.lowpassA is not None:
            im=lowpass(im,effective_pixel,args.lowpassA)
        base=autocontrast(im) if not args.no_autocontrast else (im-im.min())/(im.max()-im.min()+1e-8)
        # CryoSPARC thumbnail display coordinates require a vertical flip of
        # the micrograph background. Keep overlay coordinates and templates
        # unchanged, as requested.
        base=np.flipud(base).copy()
        original_canvas=Image.fromarray(np.uint8(base*255),'L').convert('RGBA')
        canvas=original_canvas.copy()
        templates_canvas=Image.new('RGBA',canvas.size,(0,0,0,0))
        subset_canvases=[Image.new('RGBA',canvas.size,(0,0,0,0)) for _ in allowed_sets]
        coloured_templates_canvas=Image.new('RGBA',canvas.size,(0,0,0,0))
        for r in prs:
            cmap=source_data[r['_source_index']][2]; stack=source_data[r['_source_index']][3]
            cid=int(r['_rlnClassNumber']); idx=cmap.get(cid)
            if not idx: continue
            c=stack[idx-1]
            # Binary mask from the original MRC intensity, then pad/dilate
            # in template pixels. This is deliberately not local matching.
            alpha_mask=c > args.mask_threshold
            labels,ncomp=label(alpha_mask, structure=np.ones((3,3),dtype=bool))
            if ncomp:
                counts=np.bincount(labels.ravel()); counts[0]=0
                alpha_mask=(labels == counts.argmax())
            if args.mask_pad > 0:
                yy,xx=np.ogrid[-args.mask_pad:args.mask_pad+1,-args.mask_pad:args.mask_pad+1]
                alpha_mask=binary_dilation(alpha_mask, structure=(xx*xx+yy*yy <= args.mask_pad**2))
            if not args.no_invert_classes: c=-c
            c=(c-c.min())/(c.max()-c.min()+1e-8)
            # class stack is 2.390625 A/px; thumbnail is .85 A/px
            size=max(1,round(c.shape[0]*2.390625/args.Apix*scale)); c=Image.fromarray(np.uint8(c*255),'L').resize((size,size),Image.Resampling.BICUBIC)
            c=c.rotate(float(r.get('_rlnAnglePsi','0')),resample=Image.Resampling.BICUBIC,expand=False)
            alpha=Image.fromarray(np.uint8(alpha_mask)*255,'L').resize((size,size),Image.Resampling.NEAREST)
            # The transparency mask must undergo exactly the same geometric
            # transform as the class image, otherwise the mask remains
            # axis-aligned while the template rotates inside it.
            alpha=alpha.rotate(float(r.get('_rlnAnglePsi','0')),resample=Image.Resampling.NEAREST,expand=False)
            alpha_arr=np.asarray(alpha,dtype=np.float32)/255.
            if args.mask_smoothing>0: alpha_arr=gaussian_filter(alpha_arr,args.mask_smoothing)
            alpha=Image.fromarray(np.uint8(np.clip(alpha_arr*255*args.opacity,0,255)),'L')
            c.putalpha(alpha)
            # Apply the inverse 2D alignment shift as well as the picked
            # coordinate. STAR origin shifts are in Angstroms; thumbnails
            # are rendered at args.pixel Angstroms/pixel.
            x=round(float(r['_x_frac'])*canvas.width-size/2)
            y=round(float(r['_y_frac'])*canvas.height-size/2)
            layer=Image.new('RGBA',canvas.size); layer.paste(c,(x,y),c.getchannel('A'))
            # Colour only the combined overlay. Preserve template shading:
            # dark class density remains dark, while the background is tinted.
            subset_id=r['_subset_indices'][0] if r['_subset_indices'] else None
            if allowed_sets and any(x is not None for x in subset_colors) and subset_id is not None and subset_colors[subset_id] is not None and subset_modes[subset_id] != 'recolour':
                # Keep the grayscale template and add a 2-pixel outer edge.
                mb=alpha_arr>0.5
                # Straddle the smoothed template edge: using only the outer
                # dilation leaves a gap between the antialiased template and
                # the coloured border.
                border=binary_dilation(mb,iterations=args.border_thickness) & ~binary_erosion(mb,iterations=args.border_thickness)
                if subset_modes[subset_id]=='border_dashed':
                    by,bx=np.indices(border.shape); ang=(np.arctan2(by-size/2,bx-size/2)+np.pi)%(2*np.pi)
                    arc=ang*size/2; border &= (arc % 40 < 24)
                ba=Image.fromarray(np.uint8(border)*255,'L')
                bc=Image.new('RGBA',(size,size),subset_colors[subset_id]+(0,)); bc.putalpha(ba)
                clayer=Image.new('RGBA',canvas.size); clayer.paste(c,(x,y),c.getchannel('A'))
                blayer=Image.new('RGBA',canvas.size); blayer.paste(bc,(x,y),ba)
                canvas=Image.alpha_composite(canvas,clayer); canvas=Image.alpha_composite(canvas,blayer)
                coloured_templates_canvas=Image.alpha_composite(coloured_templates_canvas,clayer)
                coloured_templates_canvas=Image.alpha_composite(coloured_templates_canvas,blayer)
            elif allowed_sets and any(x is not None for x in subset_colors):
                rgb=np.asarray(c.convert('L'),dtype=np.float32)/255.
                colour=subset_colors[subset_id] if subset_id is not None and subset_colors[subset_id] is not None else (128,128,128)
                tint=np.asarray(colour,dtype=np.float32)/255.
                rgb_col=np.uint8(np.clip(rgb[...,None]*(0.35+0.65*tint)*255,0,255))
                rgba=np.dstack([np.uint8(np.clip(rgb[...,None]*(0.35+0.65*tint[None,None,:])*255,0,255)),np.asarray(c.getchannel('A'))])
                cc=Image.fromarray(np.uint8(rgba),'RGBA')
                clayer=Image.new('RGBA',canvas.size); clayer.paste(cc,(x,y),cc.getchannel('A'))
                canvas=Image.alpha_composite(canvas,clayer)
                coloured_templates_canvas=Image.alpha_composite(coloured_templates_canvas,clayer)
            else:
                canvas=Image.alpha_composite(canvas,layer)
            templates_canvas=Image.alpha_composite(templates_canvas,layer)
            for si in r['_subset_indices']:
                subset_canvases[si]=Image.alpha_composite(subset_canvases[si],layer)
        stem=p.stem.replace('_thumb_@1x','')
        if not allowed_sets:
            files=[(stem+'_1_repositioned.png',canvas.convert('RGB')),
                   (stem+'_2_original.png',original_canvas.convert('RGB')),
                   (stem+'_3_templates.png',templates_canvas)]
        else:
            files=[(stem+'_1_original.png',original_canvas.convert('RGB')),
                   (stem+'_2_repositioned_all.png',canvas.convert('RGB'))]
            for i,(subset_label,sc) in enumerate(zip(subset_labels,subset_canvases),3):
                files.append((stem+f'_{i}_{subset_label}.png',sc))
            files.append((stem+f'_{len(files)+1}_all_coloured_templates.png',coloured_templates_canvas))
        # CryoSPARC thumbnail coordinates use the opposite Y origin from the
        # display convention used here. Flip the completed render so the
        # micrograph, overlays, and transparent template layer remain in the
        # same coordinate frame. Do not alter pose angles separately.
        for fn,img in files: img.transpose(Image.Transpose.FLIP_TOP_BOTTOM).save(out/fn)
        report.append((k,len(prs),str(out/files[0][0])))
    with open(out/'report.csv','w',newline='') as f:
        w=csv.writer(f); w.writerow(['micrograph_key','particles','output']); w.writerows(report)
    print(f'Wrote {len(report)} images from {sum(x[1] for x in report)} matched particles; {len(thumbs)} thumbnails available')
    print(f'Unmatched particle micrograph keys: {len(set(mic_key(r["_rlnMicrographName"]) for r in rows)-set(thumbs))}')
if __name__=='__main__': main()
