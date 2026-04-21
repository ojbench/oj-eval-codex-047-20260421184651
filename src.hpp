#ifndef PPCA_SRC_HPP
#define PPCA_SRC_HPP

#include "math.h"
#include <algorithm>

class Monitor; // Forward declaration; provided by the judge codebase

class Controller {

public:
    Controller(const Vec &_pos_tar, double _v_max, double _r, int _id, Monitor *_monitor) {
        pos_tar = _pos_tar;
        v_max = _v_max;
        r = _r;
        id = _id;
        monitor = _monitor;
    }

    void set_pos_cur(const Vec &_pos_cur) {
        pos_cur = _pos_cur;
    }

    void set_v_cur(const Vec &_v_cur) {
        v_cur = _v_cur;
    }

private:
    int id = 0;
    Vec pos_tar;
    Vec pos_cur;
    Vec v_cur;
    double v_max = 0, r = 0;
    Monitor *monitor = nullptr;

    // Compute desired velocity towards target with speed cap and no overshoot
    Vec desired_to_target() const {
        Vec to_tar = pos_tar - pos_cur;
        double dist = to_tar.norm();
        if (dist <= 1e-8) return Vec();
        double speed = std::min(v_max, dist / TIME_INTERVAL);
        return to_tar.normalize() * speed;
    }

    // Check if moving with velocity v_self is collision-free within TIME_INTERVAL
    // against other robots, assuming others keep their last known velocity.
    bool collision_free_with(const Vec &v_self) const;

    // Yield rule: give priority to smaller-id neighbors within a guard radius.
    bool should_yield() const;

public:

    Vec get_v_next() {
        // If already sufficiently close, stop to avoid jittering around the target
        if ((pos_cur - pos_tar).norm_sqr() <= 1.5 * 1.5 * 1e-4) {
            return Vec();
        }

        // Basic desired velocity towards the target
        Vec v_des = desired_to_target();

        // Simple priority-based yielding to avoid deadlocks/collisions
        if (should_yield()) {
            return Vec();
        }

        // Scale down speed via binary search until predicted collision-free
        double lo = 0.0, hi = 1.0;
        for (int it = 0; it < 25; ++it) {
            double mid = (lo + hi) * 0.5;
            Vec v_try = v_des * mid;
            if (collision_free_with(v_try)) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        Vec v_final = v_des * lo;
        // Final safety clamp to v_max
        double sp = v_final.norm();
        if (sp > v_max) v_final = v_final * (v_max / sp);
        return v_final;
    }
};


/////////////////////////////////
/// Helper implementations
/////////////////////////////////

#include "monitor.h"

inline bool Controller::collision_free_with(const Vec &v_self) const {
    int n = monitor->get_robot_number();
    for (int j = 0; j < n; ++j) {
        if (j == id) continue;
        Vec pj = monitor->get_pos_cur(j);
        Vec vj = monitor->get_v_cur(j);
        double rj = monitor->get_r(j);

        Vec delta_pos = pos_cur - pj;
        double d = delta_pos.norm();
        double delta_r = r + rj;
        if (d <= delta_r + 1e-6) {
            // Too close already — stay put to avoid worsening overlap
            return false;
        }

        Vec delta_v = v_self - vj;
        double project = delta_pos.dot(delta_v);
        if (project >= 0) {
            // Relative motion moving apart — safe against this neighbor
            continue;
        }
        project /= -delta_v.norm();
        double min_dis_sqr;
        if (project < delta_v.norm() * TIME_INTERVAL) {
            min_dis_sqr = delta_pos.norm_sqr() - project * project;
        } else {
            min_dis_sqr = (delta_pos + delta_v * TIME_INTERVAL).norm_sqr();
        }
        if (min_dis_sqr <= delta_r * delta_r - EPSILON) {
            return false;
        }
    }
    return true;
}

inline bool Controller::should_yield() const {
    int n = monitor->get_robot_number();
    // Guard radius: sum of radii + a buffer proportional to our max movement
    double guard = v_max * TIME_INTERVAL + 2.0 * r + 1e-3;
    for (int j = 0; j < n; ++j) {
        if (j == id) continue;
        double rj = monitor->get_r(j);
        Vec pj = monitor->get_pos_cur(j);
        double d = (pos_cur - pj).norm();
        if (d <= (r + rj) + guard) {
            if (j < id) return true; // yield to smaller id in proximity
        }
    }
    return false;
}

#endif // PPCA_SRC_HPP

