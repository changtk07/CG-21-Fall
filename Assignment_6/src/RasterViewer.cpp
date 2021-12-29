#include "SDLViewer.h"

#include <Eigen/Core>
#include <Eigen/Dense>

#include <functional>
#include <iostream>
#include <climits>

#include "raster.h"

// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION // Do not include this line twice in your project!
#include "stb_image_write.h"

unsigned select_triangle(const Eigen::Vector4f &p, const std::vector<VertexAttributes> &vertices, const UniformAttributes &uniform) {
  unsigned select = -1;
  for (unsigned i = 0; i < vertices.size()/3; ++i) {
    Eigen::Vector4f a = uniform.view * vertices[i*3+0].transform * vertices[i*3+0].position;
    Eigen::Vector4f b = uniform.view * vertices[i*3+1].transform * vertices[i*3+1].position;
    Eigen::Vector4f c = uniform.view * vertices[i*3+2].transform * vertices[i*3+2].position;
    Eigen::Vector4f u = b-a, v = c-a, d(0,0,1,0);
    Eigen::Matrix4f A = Eigen::Matrix4f::Zero();
    A << -u, -v, d;
    Eigen::Vector4f r = a-p;
    Eigen::Vector4f x = A.colPivHouseholderQr().solve(r);
    if (x[0] >= 0 && x[1] >= 0 && x[0]+x[1] <= 1) {
      select = i;
    }
  }
  return select;
}

unsigned closest_vertex(const Eigen::Vector4f &p, const std::vector<VertexAttributes> &vertices, const UniformAttributes &uniform) {
  unsigned cloest = 0;
  float mindist = std::numeric_limits<float>::max();
  for (unsigned i = 0; i < vertices.size(); ++i) {
    float dist = (uniform.view * vertices[i].transform * vertices[i].position - p).squaredNorm();
    if (dist < mindist) {
      mindist = dist;
      cloest = i;
    }
  }
  return cloest;
}

Eigen::Matrix4f scale(float pct) {
  Eigen::Matrix4f M;
  M <<
  1+pct, 0, 0, 0,
  0, 1+pct, 0, 0,
  0, 0, 1+pct, 0,
  0, 0, 0, 1;
  return M;
}

Eigen::Matrix4f rot2d(float d, bool radian=false) {
  if (!radian) d = d / 180 * M_PI;
  Eigen::Matrix4f M;
  M <<
  std::cos(d), -std::sin(d), 0, 0,
  std::sin(d), std::cos(d), 0, 0,
  0, 0, 1, 0,
  0, 0, 0, 1;
  return M;
}

unsigned comb(unsigned top, unsigned bottom) {
  unsigned x = 1, y = 1;
  for (unsigned i=0; i <bottom; ++i) {
    x *= top-i;
    y *= bottom-i;
  }
  return x/y;
}

int main(int argc, char *args[]) {
  int width = 800;
  int height = 800;

  std::vector<Eigen::Vector4f> colors {
    Eigen::Vector4f(1,0,0,1),
    Eigen::Vector4f(0,1,0,1),
    Eigen::Vector4f(0,0,1,1),
    Eigen::Vector4f(1,1,0,1),
    Eigen::Vector4f(1,0,1,1),
    Eigen::Vector4f(0,1,1,1),
    Eigen::Vector4f(1,0.5,0,1),
    Eigen::Vector4f(0.5,0.5,0.5,1),
    Eigen::Vector4f(0,0,0,1)
  };

  enum VIWER_MODE {
    DEFAULT,
    INSERT,
    TRANSLATE,
    DELETE,
    COLOR
  } mode, prev;

  VertexAttributes tmpv;
  unsigned select = -1;
  bool mouse_click = false;
  bool mouse_hold = false;
  float highlight = 0.3;

  std::vector<Eigen::Matrix4f> keyframes;
  std::vector<Eigen::Matrix4f> step(3);
  const unsigned K = 50;
  unsigned k = K;
  unsigned target = -1;
  unsigned currframe = -1;
  bool bezier = false;
  bool animating = false;
  
	Eigen::Matrix<FrameBufferAttributes,Eigen::Dynamic,Eigen::Dynamic> frameBuffer(width, height); // The Framebuffer storing the image rendered by the rasterizer
	UniformAttributes uniform; // Global Constants
	Program program; // Basic rasterization program

	// The vertex shader
	program.VertexShader = [&](const VertexAttributes& va, const UniformAttributes& uniform) {
    VertexAttributes out;
    out.position = uniform.view * va.transform * va.position;
    out.color = va.color.cwiseMax(0);
    out.color = va.color.cwiseMin(1);
		return out;
	};

	// The fragment shader
	program.FragmentShader = [](const VertexAttributes& va, const UniformAttributes& uniform) {
		return FragmentAttributes(va.color(0),va.color(1),va.color(2));
	};

	// The blending shader converts colors between 0 and 1 to uint8
	program.BlendingShader = [](const FragmentAttributes& fa, const FrameBufferAttributes& previous) {
		return FrameBufferAttributes(fa.color[0]*255,fa.color[1]*255,fa.color[2]*255,fa.color[3]*255);
	};

	// vertices
	std::vector<VertexAttributes> vertices;
  std::vector<VertexAttributes> dirty;

  // Initialize the viewer and the corresponding callbacks
  SDLViewer viewer;
  viewer.init("Viewer Example", width, height);

  viewer.mouse_move = [&](int x, int y, int xrel, int yrel) {
    if (animating) return;
    // std::cout << "LOG: mouse moved " << x << "," << y << " " << xrel << "," << yrel << "\n";
    Eigen::Vector4f p((float(x)/float(width)*2)-1, (float(height-1-y)/float(height)*2)-1, 0, 1);
    Eigen::Vector4f crel(float(xrel)/float(width)*2, -float(yrel)/float(height)*2, 0, 0);

    if (mode == INSERT) {
      tmpv.position = uniform.view.inverse() * p;
    } else if (mode==TRANSLATE && mouse_hold && select<vertices.size()) {
      Eigen::Matrix4f inv = uniform.view.inverse();
      vertices[select*3+0].transform.col(3) += inv*crel;
      vertices[select*3+1].transform.col(3) += inv*crel;
      vertices[select*3+2].transform.col(3) += inv*crel;
    }
    viewer.redraw_next = true;
  };

  viewer.mouse_pressed = [&](int x, int y, bool is_pressed, int button, int clicks) {
    if (animating) return;
    std::cout << "LOG: mouse " << (is_pressed ? "pressed" : "released") << " (" 
              << (button==1 ? "left" : button==2 ? "mid" : "right" ) << ") "
              << x << " " << y << " " << "\n";
    mouse_hold = is_pressed;
    
    Eigen::Vector4f p((float(x)/float(width)*2)-1, (float(height-1-y)/float(height)*2)-1, 0, 1);
    if (mode == INSERT && is_pressed) {
      mouse_click = true;
      tmpv.position = uniform.view.inverse() * p;
    } else if (mode == TRANSLATE) {
      if (!is_pressed && select<vertices.size()/3) {
        vertices[select*3+0].color.array() -= highlight;
        vertices[select*3+1].color.array() -= highlight;
        vertices[select*3+2].color.array() -= highlight;
      } else if (is_pressed) {
        select = select_triangle(p, vertices, uniform);
        if (target!=-1 && select!=target) select = -1;
        if (select < vertices.size()/3) {
          vertices[select*3+0].color.array() += highlight;
          vertices[select*3+1].color.array() += highlight;
          vertices[select*3+2].color.array() += highlight;
        }
      }
    } else if (mode == DELETE && is_pressed) {
      select = select_triangle(p, vertices, uniform);
      if (select < vertices.size()/3) {
        vertices.erase(vertices.begin()+select*3, vertices.begin()+select*3+3);
      }
    } else if (mode == COLOR && is_pressed) {
      select = closest_vertex(p, vertices, uniform);
    }
    viewer.redraw_next = true;
  };

  viewer.mouse_wheel = [&](int dx, int dy, bool is_direction_normal) {
  };

  viewer.key_pressed = [&](char key, bool is_pressed, int modifier, int repeat) {
    if (animating) return;
    std::cout << "LOG: key " << (is_pressed ? "pressed" : "released") << " " << (int)key << "(" << key << ")\n";
    if (!is_pressed) return;
    prev = mode;

    if (key == 27) { // ESC
      mode = DEFAULT;
    } else if (key == 'i') {
      mode = INSERT;
    } else if (key == 'o') {
      mode = TRANSLATE;
    } else if (key == 'p') {
      mode = DELETE;
    } else if (key == 'c') {
      mode = COLOR;
    } else if ((key=='h' || key=='j' || key=='k' || key=='l') && mode==TRANSLATE && select < vertices.size()/3) {
      Eigen::Matrix4f M = key=='h' ? rot2d(10) : key=='j' ? rot2d(-10) : key=='k' ? scale(.25) : scale(-.25);
      VertexAttributes &a = vertices[select*3+0];
      VertexAttributes &b = vertices[select*3+1];
      VertexAttributes &c = vertices[select*3+2];
      Eigen::Vector4f barycenter = (a.position+b.position+c.position)/3;
      Eigen::Matrix4f tr1 = Eigen::Matrix4f::Identity();
      Eigen::Matrix4f tr2 = Eigen::Matrix4f::Identity();
      tr1.col(3) = -barycenter;
      tr2.col(3) = barycenter;
      tr1(3,3) = 1;
      a.transform *= tr2 * M * tr1;
      b.transform *= tr2 * M * tr1;
      c.transform *= tr2 * M * tr1;
    } else if (key >= '1' && key <= '9' && mode == COLOR && select < vertices.size()) {
      vertices[select].color = colors[key-'1'];
    } else if (key == '=' && mode == DEFAULT) {
      uniform.view *= scale(.2);
    } else if (key == '-' && mode == DEFAULT) {
      uniform.view *= scale(-.2);
    } else if (key == 'w' && mode == DEFAULT) {
      uniform.view.col(3) += Eigen::Vector4f(0, -0.4, 0, 0);
    } else if (key == 'a' && mode == DEFAULT) {
      uniform.view.col(3) += Eigen::Vector4f(0.4, 0, 0, 0);
    } else if (key == 's' && mode == DEFAULT) {
      uniform.view.col(3) += Eigen::Vector4f(0, 0.4, 0, 0);
    } else if (key == 'd' && mode == DEFAULT) {
      uniform.view.col(3) += Eigen::Vector4f(-0.4, 0, 0, 0);
    } else if (key == 'z' && mode == TRANSLATE && select < vertices.size()/3) { // add a keyframe
      target = select;
      keyframes.push_back( vertices[select*3+0].transform );
      keyframes.push_back( vertices[select*3+1].transform );
      keyframes.push_back( vertices[select*3+2].transform );
    } else if (key == 'x' && mode == TRANSLATE) { // clear keyframes
      keyframes.clear();
      target = -1;
    } else if (key=='v' && mode == TRANSLATE && !keyframes.empty() && keyframes.size()/3) { // linear
      bezier = false;
      animating = true;
      currframe = -1;
      k = K;
    } else if (key=='b' && mode == TRANSLATE && !keyframes.empty() && keyframes.size()/3) { // bezier
      bezier = true;
      animating = true;
      k = 0;
    }

    if (mode != prev) {
      if (prev==TRANSLATE && select<vertices.size()/3 && mouse_hold) {
        vertices[select*3+0].color.array() -= highlight;
        vertices[select*3+1].color.array() -= highlight;
        vertices[select*3+2].color.array() -= highlight;
      }
      select = -1;
      target = -1;
      dirty.clear();
      keyframes.clear();
    }
    viewer.redraw_next = true;
  };

  viewer.redraw = [&](SDLViewer &viewer) {
    // std::cout << "LOG: re-draw\n";
    // Clear the framebuffer
    for (unsigned i=0; i<frameBuffer.rows(); i++)
      for (unsigned j=0; j<frameBuffer.cols(); j++)
        frameBuffer(i,j).color << 255,255,255,255;
    
    std::size_t fixed = dirty.size();
    if (dirty.size() < 2) dirty.push_back(tmpv);
    else dirty.insert(dirty.end(), {dirty[1], tmpv, tmpv, dirty[0]});

    if (mouse_click && dirty.size()==6) {
      vertices.push_back(dirty[0]);
      vertices.back().color = colors[0];
      vertices.push_back(dirty[2]);
      vertices.back().color = colors[0];
      vertices.push_back(dirty[4]);
      vertices.back().color = colors[0];
      dirty.clear();
    }

    if (animating && !bezier) {
      if (k == K) {
        k = 0;
        ++currframe;
        unsigned nextframe = currframe+1;
        step[0] = (keyframes[nextframe*3+0]-keyframes[currframe*3+0]).array() / K;
        step[1] = (keyframes[nextframe*3+1]-keyframes[currframe*3+1]).array() / K;
        step[2] = (keyframes[nextframe*3+2]-keyframes[currframe*3+2]).array() / K;
      }

      vertices[target*3+0].transform = keyframes[currframe*3+0] + k*step[0];
      vertices[target*3+1].transform = keyframes[currframe*3+1] + k*step[1];
      vertices[target*3+2].transform = keyframes[currframe*3+2] + k*step[2];

      ++k;
      animating = currframe+1 < keyframes.size()/3;
      viewer.redraw_next = true;
    } else if (animating && bezier) {
      auto n = keyframes.size()/3;
      float t = float(k)/(K*(n-1));
      
      vertices[target*3+0].transform.setZero();
      vertices[target*3+1].transform.setZero();
      vertices[target*3+2].transform.setZero();

      for (unsigned i=0; i<n; ++i) {
        float coeff = comb(n-1, i) * std::pow(t, i) * std::pow(1-t, n-i-1);
        vertices[target*3+0].transform += coeff * keyframes[i*3+0];
        vertices[target*3+1].transform += coeff * keyframes[i*3+1];
        vertices[target*3+2].transform += coeff * keyframes[i*3+2];
      }
      ++k;
      animating = k <= K*(n-1);
      viewer.redraw_next = animating;
    }

    for (unsigned i=0; i<vertices.size()/3; ++i) {
      rasterize_triangles(program, uniform, {vertices[i*3+0], vertices[i*3+1], vertices[i*3+2]}, frameBuffer);
      std::vector<VertexAttributes> lines = {vertices[i*3+0], vertices[i*3+1], vertices[i*3+1], vertices[i*3+2], vertices[i*3+2], vertices[i*3+0]};
      std::for_each(lines.begin(), lines.end(), [&](VertexAttributes &v) { v.color = colors.back(); });
      rasterize_lines(program, uniform, lines, 0.6, frameBuffer);
    }

    rasterize_lines(program, uniform, dirty, 0.6, frameBuffer);

    if (!mouse_click) dirty.resize(fixed);
    mouse_click = false;

    // Buffer for exchanging data between rasterizer and sdl viewer
    Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic> R(width, height);
    Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic> G(width, height);
    Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic> B(width, height);
    Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic> A(width, height);

    for (unsigned i=0; i<frameBuffer.rows();i++) {
      for (unsigned j=0; j<frameBuffer.cols();j++) {
        R(i,frameBuffer.cols()-1-j) = frameBuffer(i,j).color(0);
        G(i,frameBuffer.cols()-1-j) = frameBuffer(i,j).color(1);
        B(i,frameBuffer.cols()-1-j) = frameBuffer(i,j).color(2);
        A(i,frameBuffer.cols()-1-j) = frameBuffer(i,j).color(3);
      }
    }
    viewer.draw_image(R, G, B, A);
  };

  viewer.launch();

  return 0;
}
