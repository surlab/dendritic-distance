# Annotation

Additional details about the required annotations. For detailed instructions on how to generate these annotations using Fiji (ImageJ) please see the PDF document "ROI selection."

## dendrite annotation file
1 file whose name includes the string "dend." that contains annotations of the dendrites and branch points

Currently this file should be produced using ImageJ to generate ROIs (if necessary it can easily be adapted to other formats). 
The order of ROIs does not matter for this file since they are well identified by ROI type
Each dendritic segment should be marked with a segmented line ROI in imageJ
Each branch point should be marked with an Oval ROI in imageJ

It is ~~recomended to~~ currently required that 3 dendrite segments connect to each Oval marked branch point, even if one of the segments contains no spines or is very tiny before moving out of focus or off screen. This could be remedied in the future by applying a threshold distance at which a segment cannot be connected to a branch. 

It is also recomended to annotate all segments in view even if they have no spines, aren't on the branch of interest, or may be axons. This will help assess accuracy if we move to an automated algorithm for identifying segments in the future. 

Please see "demo_data" for examples of "*dend.zip and "*dend.roi" dendrite annotation files



## spine annotation file
1 file whose name includes the string ".zip" and not the string ".dend" that contains the annotations of rois for each spine head and corresponding nearest dendritic segment

Currently this file should be produced using ImageJ to generate ROIs (if necessary it can easily be adapted to other formats). 
The order of ROIs in this file matters since otherwise what they correspond to is ambiguous. The order should be 
1. Spine Head #1
2. Background #1 (same shape and size as spine head, this ROI is not used by this code)
3. Dendrite #1 (The center of this ROI is used to identify which of the dendritic segments (see above) the spine should be connected to)
4. Spine head #2
5. Background #2
6. etc.

Different types of ROIs should be fine for each of these. The only recomendation is that for spines close to a branch point, do not place the "dendrite" roi so that it spans the branch point. This can lead to ambiguity and mis-attricution to the wrong dendritic segment. Place the "dendrite" ROI so that it resides entirely on the same dendritic segment that you believe the spin is connected to. 

Please see "demo_data" for examples of "*.zip" spine annotation files
