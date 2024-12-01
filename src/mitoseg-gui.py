#!/usr/bin/python3

import glob
import multiprocessing
import subprocess
import tkinter as tk
from threading import Thread
from tkinter import ttk, filedialog, messagebox

process = None


def run_binary():
    global process

    if process and process.poll() is None:
        response = messagebox.askyesno(
            'Process Running',
            'A process is already running. Do you want to terminate it and start a new one?'
        )

        if response:
            process.terminate()
            process.wait()
        else:
            return

    command = [
        './mitoseg', '--zrange',
        str(zrange_start.get()),
        str(zrange_end.get()), '--psize',
        str(psize.get())
    ]

    if pattern.get() != '':
        command.extend(['--pattern', pattern.get()])

    if roi_mode.get() == 'manual':
        command.extend([
            '--roi',
            str(roi_x.get()),
            str(roi_y.get()),
            str(roi_w.get()),
            str(roi_h.get())
        ])

    if src.get():
        command.extend(['--src', src.get()])
    if dst.get():
        command.extend(['--dst', dst.get()])
    if phase.get() != 'All':
        command.extend(['--phase', phase.get()])
    if valid.get():
        command.extend(['--valid', valid.get()])
    if thick_full.get():
        command.extend(['--thick', 'full'])
    elif thick.get():
        command.extend(['--thick', thick.get()])
    if cores.get():
        command.extend(['--cores', str(cores.get())])
    if settings_file.get() and settings_file.get() != 'None':
        command.extend(['--settings-file', settings_file.get()])

    text_box.delete('1.0', tk.END)

    process = subprocess.Popen(command,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               text=True)

    def read_output():
        for line in process.stdout:
            text_box.insert(tk.END, line)
            text_box.see(tk.END)
        for line in process.stderr:
            text_box.insert(tk.END, line)
            text_box.see(tk.END)
        process.wait()

        if process.returncode == 0 and dst.get():
            subprocess.run(['xdg-open', dst.get()])

    Thread(target=read_output, daemon=True).start()


def toggle_roi(enable):
    state = 'normal' if enable else 'disabled'
    roi_x_entry.config(state=state)
    roi_y_entry.config(state=state)
    roi_w_entry.config(state=state)
    roi_h_entry.config(state=state)


def browse_src():
    folder_path = filedialog.askdirectory()
    if folder_path:
        src.set(folder_path)


def browse_dst():
    folder_path = filedialog.askdirectory()
    if folder_path:
        dst.set(folder_path)


def toggle_thick_entry():
    if thick_full.get():
        thick_entry.config(state='disabled')
    else:
        thick_entry.config(state='normal')


root = tk.Tk()
root.title('MitoSeg v1.0')

left_frame = tk.Frame(root)
left_frame.grid(row=0, column=0, padx=10, pady=10, sticky='nsew')

right_frame = tk.Frame(root)
right_frame.grid(row=0, column=1, padx=10, pady=10, sticky='nsew')

mandatory_frame = tk.LabelFrame(left_frame, text='Mandatory Settings')
mandatory_frame.pack(fill='x', padx=5, pady=5)
mandatory_frame.grid_columnconfigure(1, weight=1)

zrange_start_label = tk.Label(mandatory_frame, text='Z-range start:')
zrange_start_label.grid(row=0, column=0, padx=5, pady=5, sticky='w')
zrange_start = tk.IntVar(value=35)
zrange_start_entry = tk.Entry(mandatory_frame, textvariable=zrange_start)
zrange_start_entry.grid(row=0, column=1, padx=5, pady=5, sticky='ew')

zrange_end_label = tk.Label(mandatory_frame, text='Z-range end:')
zrange_end_label.grid(row=1, column=0, padx=5, pady=5, sticky='w')
zrange_end = tk.IntVar(value=74)
zrange_end_entry = tk.Entry(mandatory_frame, textvariable=zrange_end)
zrange_end_entry.grid(row=1, column=1, padx=5, pady=5, sticky='ew')

psize_label = tk.Label(mandatory_frame, text='Pixel size (px/nm):')
psize_label.grid(row=2, column=0, padx=5, pady=5, sticky='w')
psize = tk.DoubleVar(value=2.2)
psize_entry = tk.Entry(mandatory_frame, textvariable=psize)
psize_entry.grid(row=2, column=1, padx=5, pady=5, sticky='ew')

pattern_label = tk.Label(mandatory_frame, text='Filename pattern:')
pattern_label.grid(row=3, column=0, padx=5, pady=5, sticky='w')
pattern = tk.StringVar(value='gap18_sub%04d.png')
pattern_entry = tk.Entry(mandatory_frame, textvariable=pattern)
pattern_entry.grid(row=3, column=1, padx=5, pady=5, sticky='ew')

roi_frame = tk.LabelFrame(left_frame, text='Region of Interest')
roi_frame.pack(fill='x', padx=5, pady=5)

roi_mode = tk.StringVar(value='auto')

auto_radio = tk.Radiobutton(roi_frame,
                            text='Auto',
                            variable=roi_mode,
                            value='auto',
                            command=lambda: toggle_roi(False))
auto_radio.grid(row=0, column=0, padx=5, pady=5, sticky='w')

manual_radio = tk.Radiobutton(roi_frame,
                              text='Manual',
                              variable=roi_mode,
                              value='manual',
                              command=lambda: toggle_roi(True))
manual_radio.grid(row=0, column=1, padx=5, pady=5, sticky='w')

roi_x_label = tk.Label(roi_frame, text='X:')
roi_x_label.grid(row=1, column=0, padx=5, pady=5, sticky='w')
roi_x = tk.StringVar()
roi_x_entry = tk.Entry(roi_frame, textvariable=roi_x, state='disabled')
roi_x_entry.grid(row=1, column=1, padx=5, pady=5, sticky='ew')

roi_y_label = tk.Label(roi_frame, text='Y:')
roi_y_label.grid(row=2, column=0, padx=5, pady=5, sticky='w')
roi_y = tk.StringVar()
roi_y_entry = tk.Entry(roi_frame, textvariable=roi_y, state='disabled')
roi_y_entry.grid(row=2, column=1, padx=5, pady=5, sticky='ew')

roi_w_label = tk.Label(roi_frame, text='Width:')
roi_w_label.grid(row=3, column=0, padx=5, pady=5, sticky='w')
roi_w = tk.StringVar()
roi_w_entry = tk.Entry(roi_frame, textvariable=roi_w, state='disabled')
roi_w_entry.grid(row=3, column=1, padx=5, pady=5, sticky='ew')

roi_h_label = tk.Label(roi_frame, text='Height:')
roi_h_label.grid(row=4, column=0, padx=5, pady=5, sticky='w')
roi_h = tk.StringVar()
roi_h_entry = tk.Entry(roi_frame, textvariable=roi_h, state='disabled')
roi_h_entry.grid(row=4, column=1, padx=5, pady=5, sticky='ew')

roi_frame.grid_columnconfigure(1, weight=1)

optional_frame = tk.LabelFrame(right_frame, text='Optional Settings')
optional_frame.pack(fill='x', padx=5, pady=5)
optional_frame.grid_columnconfigure(1, weight=1)

src_label = tk.Label(optional_frame, text='Source path:')
src_label.grid(row=0, column=0, padx=5, pady=5, sticky='w')
src = tk.StringVar()
src_entry = tk.Entry(optional_frame, textvariable=src)
src_entry.grid(row=0, column=1, padx=5, pady=5, sticky='ew')
src_button = tk.Button(optional_frame, text='Browse...', command=browse_src)
src_button.grid(row=0, column=2, padx=5, pady=5, sticky='ew')

dst_label = tk.Label(optional_frame, text='Destination path:')
dst_label.grid(row=1, column=0, padx=5, pady=5, sticky='w')
dst = tk.StringVar()
dst_entry = tk.Entry(optional_frame, textvariable=dst)
dst_entry.grid(row=1, column=1, padx=5, pady=5, sticky='ew')
dst_button = tk.Button(optional_frame, text='Browse...', command=browse_dst)
dst_button.grid(row=1, column=2, padx=5, pady=5, sticky='ew')

valid_label = tk.Label(optional_frame, text='Validity threshold:')
valid_label.grid(row=2, column=0, padx=5, pady=5, sticky='w')
valid = tk.StringVar()
valid_entry = tk.Entry(optional_frame, textvariable=valid)
valid_entry.grid(row=2, column=1, columnspan=2, padx=5, pady=5, sticky='ew')

thick_label = tk.Label(optional_frame, text='Snake z-thickness:')
thick_label.grid(row=3, column=0, padx=5, pady=5, sticky='w')
thick = tk.StringVar()
thick_entry = tk.Entry(optional_frame, textvariable=thick)
thick_entry.grid(row=3, column=1, padx=5, pady=5, sticky='ew')
thick_full = tk.BooleanVar(value=False)
thick_full_check = tk.Checkbutton(optional_frame,
                                  text='Full z-range',
                                  variable=thick_full,
                                  command=toggle_thick_entry)
thick_full_check.grid(row=3, column=2, padx=5, pady=5, sticky='w')

phase_label = tk.Label(optional_frame, text='Phase:')
phase_label.grid(row=4, column=0, padx=5, pady=5, sticky='w')
phase = tk.StringVar(value='All')
phase_entry = ttk.Combobox(optional_frame,
                           textvariable=phase,
                           values=['All', 1, 2, 3],
                           state='readonly')
phase_entry.grid(row=4, column=1, columnspan=2, padx=5, pady=5, sticky='ew')

cores_label = tk.Label(optional_frame, text='CPU cores:')
cores_label.grid(row=5, column=0, padx=5, pady=5, sticky='w')
coresList = list(range(1, multiprocessing.cpu_count() + 1))
cores = tk.IntVar(value=int(coresList[-1]))
cores_entry = ttk.Combobox(optional_frame,
                           textvariable=cores,
                           values=coresList,
                           state='readonly')
cores_entry.grid(row=5, column=1, columnspan=2, padx=5, pady=5, sticky='ew')

settings_file_label = tk.Label(optional_frame, text='Settings file:')
settings_file_label.grid(row=6, column=0, padx=5, pady=5, sticky='w')
settings_file = tk.StringVar(value='None')
settings_files = ['None'] + sorted(glob.glob('../*.yaml'))
settings_file_entry = ttk.Combobox(optional_frame,
                                   textvariable=settings_file,
                                   values=settings_files,
                                   state='readonly')
settings_file_entry.grid(row=6,
                         column=1,
                         columnspan=2,
                         padx=5,
                         pady=5,
                         sticky='ew')

output_frame = tk.LabelFrame(root, text='Output')
output_frame.grid(row=1,
                  column=0,
                  columnspan=2,
                  padx=10,
                  pady=10,
                  sticky='nsew')
text_box = tk.Text(output_frame, wrap=tk.WORD, height=10, width=80)
text_box.pack(fill='both', expand=True, padx=5, pady=5)

root.grid_rowconfigure(1, weight=1)
root.grid_columnconfigure(0, weight=1)
root.grid_columnconfigure(1, weight=1)

run_button = tk.Button(root, text='Start', command=run_binary)
run_button.grid(row=2, column=0, columnspan=2, pady=10)

root.mainloop()
