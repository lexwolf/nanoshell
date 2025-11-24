#!/bin/bash

flare=false
radius_ratio=""

show_usage() {
    echo "Usage: bash $0 [-f|--flare] <rho>"
    echo "  -f, --flare  Overlay the flare effect from the original sketch"
}

while [ $# -gt 0 ]; do
    case "$1" in
        -f|--flare)
            flare=true
            shift
            ;;
        -h|--help)
            show_usage
            exit 0
            ;;
        -*)
            echo "Unknown option: $1" >&2
            show_usage
            exit 1
            ;;
        *)
            radius_ratio="$1"
            shift
            break
            ;;
    esac
done

if [ -z "$radius_ratio" ]; then
    show_usage
    exit 1
fi

g++ ../src/eV2ex.cxx -o ../bin/eV2ex
g++ -Wall -I/usr/include/ -L/usr/local/lib ../src/rho2ome_sp.cxx -o ../bin/rho2ome_sp -lgsl -lgslcblas -lm -larmadillo
resu=($(./../bin/rho2ome_sp "$radius_ratio"))

inner_circle_color=$(./../bin/eV2ex "${resu[0]}")
echo "${resu[0]} $inner_circle_color"

# Set the colors for the circles (in hex format)
outer_circle_color="#C0C0C0CC"    # Silver
contour_color="#000000"           # Black
light_color="#FFFFFF"             # White

# Set the image dimensions
image_width=600
image_height=500

# Calculate the necessary canvas size to accommodate the circles
canvas_width=$((image_width > image_height ? image_width : image_height))
canvas_height=$canvas_width

# Calculate the radii based on the given ratio
outer_radius=$(awk "BEGIN{ print $canvas_width / 2 }")
inner_radius=$(awk "BEGIN{ print $outer_radius * $radius_ratio }")

# Calculate the offset for centering the circles on the canvas
offset_x=$((canvas_width / 2))
offset_y=$((canvas_height / 2))

# Define the light radius and sigma for the first spot
light_radius=5
light_sigma=50

# Define the parameters for the second spot
second_spot_radius=$((outer_radius / 6))
second_spot_blur_radius=$((second_spot_radius / 2))

second_spot_x=$((canvas_width  - 4*second_spot_radius))
second_spot_y=$((canvas_height - 3*second_spot_radius))

# Generate the outer circle
convert -size "${canvas_width}x${canvas_height}" xc:none \
  -fill "$outer_circle_color" \
  -stroke "$contour_color" -strokewidth 3 \
  -draw "translate $offset_x,$offset_y circle 0,0 0,$outer_radius" \
  \( +clone -background "$light_color" -blur 0x80 -shadow ${light_radius}x${light_sigma}+0+0 +repage \) \
  -composite \
  outer_circle.png

# Generate the inner circle
convert -size "${canvas_width}x${canvas_height}" xc:none \
  -fill "$inner_circle_color" \
  -stroke "$contour_color" -strokewidth 3 \
  -draw "translate $offset_x,$offset_y circle 0,0 0,$inner_radius" \
  inner_circle.png

# Generate the third circle (same size as the outer circle, no transparency yet)
convert -size "${canvas_width}x${canvas_height}" xc:none \
  -fill "$inner_circle_color" \
  -stroke "$contour_color" -strokewidth 3 \
  -draw "translate $offset_x,$offset_y circle 0,0 0,$outer_radius" \
  third_circle.png

# Apply transparency to the third circle (30% opacity for homogeneous look)
convert third_circle.png -alpha set -background none -channel A -evaluate multiply 0.3 +channel third_circle_transparent.png

# Generate the second spot as a solid white circle with blurred borders
convert -size "${canvas_width}x${canvas_height}" xc:none \
  -fill "$light_color" \
  -draw "translate $second_spot_x,$second_spot_y circle 0,0 $second_spot_radius,0" \
  -blur 0x$second_spot_blur_radius \
  second_spot.png

# Composite the circles and the light effects into a single image
convert outer_circle.png third_circle_transparent.png -gravity center -composite output1.png
convert output1.png inner_circle.png -gravity center -composite output2.png
convert output2.png second_spot.png -gravity center -composite nanoshell_base.png

final_output="../img/nanoshell.png"
if [ "$flare" = true ]; then
  composite -gravity center ../img/input/flare.png nanoshell_base.png "$final_output"
else
  mv nanoshell_base.png "$final_output"
fi

# Cleanup temporary files
rm -f outer_circle.png inner_circle.png third_circle.png third_circle_transparent.png output1.png output2.png second_spot.png
rm -f nanoshell_base.png

echo "Image generated successfully!"
