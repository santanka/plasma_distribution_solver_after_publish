from PIL import Image
from pdf2image import convert_from_path
import os

def crop_margins(img, threshold=200):
    # Convert image to grayscale
    img_gray = img.convert("L")
    
    # Convert image to binary using threshold
    binary = img_gray.point(lambda p: p > threshold and 255)
    
    # Get image data
    binary_data = list(binary.getdata())
    width, height = img.size
    
    # Find top margin
    for y in range(height):
        for x in range(width):
            if binary_data[y * width + x] == 0:
                top_margin = y
                break
        else:
            continue
        break
    
    # Find bottom margin
    for y in range(height - 1, -1, -1):
        for x in range(width):
            if binary_data[y * width + x] == 0:
                bottom_margin = y
                break
        else:
            continue
        break
    
    # Find left margin
    for x in range(width):
        for y in range(height):
            if binary_data[y * width + x] == 0:
                left_margin = x
                break
        else:
            continue
        break
    
    # Find right margin
    for x in range(width - 1, -1, -1):
        for y in range(height):
            if binary_data[y * width + x] == 0:
                right_margin = x
                break
        else:
            continue
        break
    
    # Crop the image
    img_cropped = img.crop((left_margin, top_margin, right_margin + 1, bottom_margin + 1))
    
    return img_cropped

def process_image(input_path, output_path):
    img = Image.open(input_path)
    cropped_img = crop_margins(img)
    cropped_img.save(output_path)

def process_pdf(input_path, output_path):
    # Convert PDF pages to images
    pages = convert_from_path(input_path)
    
    cropped_pages = []
    for page in pages:
        cropped_page = crop_margins(page)
        cropped_pages.append(cropped_page)
    
    # Save cropped images back to a new PDF
    cropped_pages[0].save(output_path, save_all=True, append_images=cropped_pages[1:])


def main(input_path, output_path):
    _, ext = os.path.splitext(input_path)

    if ext.lower() == ".png":
        process_image(input_path, output_path)
    elif ext.lower() == ".pdf":
        process_pdf(input_path, output_path)
    else:
        print("Unsupported file type:", ext)


grid_ionosphere_middle = 14
grid_middle_magnetosphere = 109
grid_fix = 175

BC_number = 16
min_number = 119

dir_name = f'/mnt/j/plasma_distribution_solver_after_publish/Earth_L_10_Imajo/alpha_perp_12_parallel_12/grid_{str(grid_ionosphere_middle).zfill(3)}_{str(grid_middle_magnetosphere).zfill(3)}_{str(grid_fix).zfill(3)}/'
dir_BC_name = f'boundary_condition_{str(BC_number)}/'

path_dir = f'{dir_name}{dir_BC_name}plot/'

png_file = f'{path_dir}numberdensity_electrostaticpotential_pressure_BC{str(BC_number)}_min{str(min_number).zfill(3)}.png'
pdf_file = f'{path_dir}numberdensity_electrostaticpotential_pressure_BC{str(BC_number)}_min{str(min_number).zfill(3)}.pdf'

main(png_file, png_file)
main(pdf_file, pdf_file)