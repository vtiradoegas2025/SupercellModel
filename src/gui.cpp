#ifdef ENABLE_GUI
#include <SFML/Graphics.hpp>
#include <vector>
#include <cmath>
#include <iostream>
#include <cstdio>
#include <algorithm>
#ifdef EXPORT_NPY
#include <fstream>
#include <string>
#endif
#include "simulation.hpp"

// External field dimensions come from simulation.hpp
const double CELL_SIZE = 4.0;  // pixels per grid cell

// === Tracer field linkage ===
// In integrated build, tracer is defined in equations.cpp via simulation.hpp externs.

// step_dynamics is declared in simulation.hpp

// Normalize for rendering [0, 255]
sf::Color tracerColor(double val) 
{
    val = std::clamp(val, 0.0, 1.0);
    int intensity = static_cast<int>(val * 255);
    return sf::Color(intensity, intensity, 255);  // Blue hue tracer
}

// === GUI Render Loop ===
void run_gui(int autoExportMs)
{
    int width = NR * CELL_SIZE;
    int height = NZ * CELL_SIZE;
    const int hudHeight = 48; // pixels for control bar at bottom

#ifdef EXPORT_NPY
    // Helper to write current theta slice to a NumPy .npy file (float32, shape (NZ, NR))
    auto save_theta_slice_npy = [&](int thetaIndex, const std::string& filename)
    {
        std::vector<float> buf;
        buf.resize(static_cast<size_t>(NR) * static_cast<size_t>(NZ));
        size_t idx = 0;
        for (int k = 0; k < NZ; ++k)
        {
            for (int i = 0; i < NR; ++i)
            {
                double v = std::clamp((double)tracer[i][thetaIndex][k], 0.0, 1.0);
                buf[idx++] = static_cast<float>(v);
            }
        }

        // Build .npy header v1.0
        std::string header_dict = "{'descr': '<f4', 'fortran_order': False, 'shape': (" + std::to_string(NZ) + ", " + std::to_string(NR) + "), }";
        size_t header_len = header_dict.size() + 1; // newline
        const size_t preamble = 6 + 2 + 2; // magic + ver + header_len
        size_t total = preamble + header_len;
        size_t padding = (16 - (total % 16)) % 16;
        header_len += padding;

        std::ofstream out(filename, std::ios::binary);
        if (!out)
        {
            std::cerr << "Failed to open file for writing: " << filename << "\n";
            return;
        }
        // magic
        out.write("\x93NUMPY", 6);
        out.put(static_cast<char>(1));
        out.put(static_cast<char>(0));
        // header length (little-endian uint16)
        uint16_t hl = static_cast<uint16_t>(header_len);
        char lenb[2];
        lenb[0] = static_cast<char>(hl & 0xFF);
        lenb[1] = static_cast<char>((hl >> 8) & 0xFF);
        out.write(lenb, 2);
        // header dict + padding + newline
        out.write(header_dict.c_str(), static_cast<std::streamsize>(header_dict.size()));
        for (size_t i = 0; i < header_len - (header_dict.size() + 1); ++i) out.put(' ');
        out.put('\n');
        // data
        out.write(reinterpret_cast<const char*>(buf.data()), static_cast<std::streamsize>(buf.size() * sizeof(float)));
        out.close();
        std::cout << "Saved slice to " << filename << " (shape=(" << NZ << "," << NR << "))\n";
    };
#endif

    sf::RenderWindow window(sf::VideoMode(sf::Vector2u(width, height + hudHeight)), "Tornado Tracer Field");
    window.setFramerateLimit(60);

    int thetaIndex = 0; // current theta slice to render
    bool paused = false;
    const int speedLevels[] = {1, 5, 20};
    int speedIndex = 0; // 0->1x, 1->5x, 2->20x
    int stepsPerFrame = speedLevels[speedIndex];
    const int gridSpacingCells = 10; // draw grid every N cells

    while (window.isOpen()) 
    {
        while (auto event = window.pollEvent()) 
        {
            if (event->is<sf::Event::Closed>())
                window.close();

            if (const auto* key = event->getIf<sf::Event::KeyPressed>())
            {
                if (key->scancode == sf::Keyboard::Scancode::Left)
                {
                    thetaIndex = (thetaIndex - 1 + NTH) % NTH;
                }
                else if (key->scancode == sf::Keyboard::Scancode::Right)
                {
                    thetaIndex = (thetaIndex + 1) % NTH;
                }
                else if (key->scancode == sf::Keyboard::Scancode::Space)
                {
                    paused = !paused;
                }
                else if (key->scancode == sf::Keyboard::Scancode::R)
                {
                    initialize();
                }
                else if (key->scancode == sf::Keyboard::Scancode::F)
                {
                    speedIndex = (speedIndex + 1) % 3;
                    stepsPerFrame = speedLevels[speedIndex];
                }
#ifdef EXPORT_NPY
                else if (key->scancode == sf::Keyboard::Scancode::S)
                {
                    // Save to data/ if present, else project root
                    std::string filename = "data/tracer_slice_th" + std::to_string(thetaIndex) + ".npy";
                    save_theta_slice_npy(thetaIndex, filename);
                }
#endif
            }

            if (const auto* mouse = event->getIf<sf::Event::MouseButtonPressed>())
            {
                if (mouse->button == sf::Mouse::Button::Left)
                {
                    float mx = static_cast<float>(mouse->position.x);
                    float my = static_cast<float>(mouse->position.y);
                    // Button layout
                    float barY = static_cast<float>(height);
                    float x = 10.0f;
                    float y = barY + 6.0f;
                    float btnH = static_cast<float>(hudHeight - 12);
                    float btnW = btnH; // square buttons

                    sf::FloatRect rewindBtn(sf::Vector2f(x, y), sf::Vector2f(btnW, btnH));
                    sf::FloatRect playPauseBtn(sf::Vector2f(x + (btnW + 10.0f), y), sf::Vector2f(btnW, btnH));
                    sf::FloatRect ffBtn(sf::Vector2f(x + 2.0f * (btnW + 10.0f), y), sf::Vector2f(btnW, btnH));

                    if (rewindBtn.contains(sf::Vector2f(mx, my)))
                    {
                        initialize(); // reset to initial conditions
                    }
                    else if (playPauseBtn.contains(sf::Vector2f(mx, my)))
                    {
                        paused = !paused;
                    }
                    else if (ffBtn.contains(sf::Vector2f(mx, my)))
                    {
                        speedIndex = (speedIndex + 1) % 3;
                        stepsPerFrame = speedLevels[speedIndex];
                    }
                }
            }
        }

        window.clear();

        // Advance simulation (integrated build)
        if (!paused)
        {
            for (int s = 0; s < stepsPerFrame; ++s)
            {
                step_dynamics();
            }
        }

#ifdef EXPORT_NPY
        // Auto-export current theta slice periodically if enabled
        static sf::Clock exportClock;
        if (autoExportMs > 0)
        {
            if (exportClock.getElapsedTime().asMilliseconds() >= autoExportMs)
            {
                exportClock.restart();
                std::string tmp = "data/.tracer_slice_th" + std::to_string(thetaIndex) + ".npy.tmp";
                std::string fin = "data/tracer_slice_th" + std::to_string(thetaIndex) + ".npy";
                save_theta_slice_npy(thetaIndex, tmp);
                std::rename(tmp.c_str(), fin.c_str());
            }
        }
#endif

        // Render tracer field as pixels
        for (int i = 0; i < NR; ++i) 
        {
            for (int k = 0; k < NZ; ++k) 
            {
                sf::RectangleShape cell(sf::Vector2f(CELL_SIZE, CELL_SIZE)); 
                cell.setFillColor(tracerColor(tracer[i][thetaIndex][k]));
                cell.setPosition(sf::Vector2f(i * CELL_SIZE, (height - (k + 1) * CELL_SIZE))); // flip z-axis
                window.draw(cell);
            }
        }

        // Grid overlay
        sf::Color gridColor(255, 255, 255, 40);
        // Vertical lines
        for (int i = 0; i <= NR; i += gridSpacingCells)
        {
            sf::RectangleShape line(sf::Vector2f(1.0f, static_cast<float>(height)));
            line.setPosition(sf::Vector2f(i * CELL_SIZE, 0.0f));
            line.setFillColor(gridColor);
            window.draw(line);
        }
        // Horizontal lines
        for (int k = 0; k <= NZ; k += gridSpacingCells)
        {
            sf::RectangleShape line(sf::Vector2f(static_cast<float>(width), 1.0f));
            line.setPosition(sf::Vector2f(0.0f, height - k * CELL_SIZE));
            line.setFillColor(gridColor);
            window.draw(line);
        }

        // HUD bar background
        sf::RectangleShape hud(sf::Vector2f(static_cast<float>(width), static_cast<float>(hudHeight)));
        hud.setPosition(sf::Vector2f(0.0f, static_cast<float>(height)));
        hud.setFillColor(sf::Color(20, 20, 20, 220));
        window.draw(hud);

        // Buttons
        float x = 10.0f;
        float y = static_cast<float>(height) + 6.0f;
        float btnH = static_cast<float>(hudHeight - 12);
        float btnW = btnH; // square buttons

        // Rewind (double left triangle)
        {
            sf::ConvexShape tri1(3);
            tri1.setPoint(0, sf::Vector2f(x + btnW, y));
            tri1.setPoint(1, sf::Vector2f(x, y + btnH / 2));
            tri1.setPoint(2, sf::Vector2f(x + btnW, y + btnH));
            tri1.setFillColor(sf::Color(200, 200, 200));
            window.draw(tri1);

            sf::ConvexShape tri2(3);
            float x2 = x + btnW * 0.6f;
            tri2.setPoint(0, sf::Vector2f(x2 + btnW, y));
            tri2.setPoint(1, sf::Vector2f(x2, y + btnH / 2));
            tri2.setPoint(2, sf::Vector2f(x2 + btnW, y + btnH));
            tri2.setFillColor(sf::Color(200, 200, 200));
            window.draw(tri2);
        }

        // Play/Pause
        {
            float bx = x + (btnW + 10.0f);
            if (paused)
            {
                // Play icon (triangle)
                sf::ConvexShape play(3);
                play.setPoint(0, sf::Vector2f(bx, y));
                play.setPoint(1, sf::Vector2f(bx, y + btnH));
                play.setPoint(2, sf::Vector2f(bx + btnW, y + btnH / 2));
                play.setFillColor(sf::Color(0, 200, 120));
                window.draw(play);
            }
            else
            {
                // Pause icon (two bars)
                sf::RectangleShape bar1(sf::Vector2f(btnW * 0.35f, btnH));
                bar1.setPosition(sf::Vector2f(bx, y));
                bar1.setFillColor(sf::Color(0, 200, 120));
                window.draw(bar1);
                sf::RectangleShape bar2(sf::Vector2f(btnW * 0.35f, btnH));
                bar2.setPosition(sf::Vector2f(bx + btnW * 0.6f, y));
                bar2.setFillColor(sf::Color(0, 200, 120));
                window.draw(bar2);
            }
        }

        // Fast forward (double right triangle), with color indicating speed level
        {
            float bx = x + 2.0f * (btnW + 10.0f);
            sf::Color c = (speedIndex == 0) ? sf::Color(200, 200, 200)
                           : (speedIndex == 1) ? sf::Color(255, 200, 120)
                                               : sf::Color(255, 120, 120);

            sf::ConvexShape tri1(3);
            tri1.setPoint(0, sf::Vector2f(bx, y));
            tri1.setPoint(1, sf::Vector2f(bx, y + btnH));
            tri1.setPoint(2, sf::Vector2f(bx + btnW, y + btnH / 2));
            tri1.setFillColor(c);
            window.draw(tri1);

            sf::ConvexShape tri2(3);
            float bx2 = bx + btnW * 0.4f;
            tri2.setPoint(0, sf::Vector2f(bx2, y));
            tri2.setPoint(1, sf::Vector2f(bx2, y + btnH));
            tri2.setPoint(2, sf::Vector2f(bx2 + btnW, y + btnH / 2));
            tri2.setFillColor(c);
            window.draw(tri2);
        }

        window.display();
    }
}

#endif // ENABLE_GUI